#!/usr/bin/env python

### Created by Bruno de Medeiros (souzademedeiros@fas.harvard.edu), starting on 19-oct-2016
### The purpose of the script is to clean 18S and 28S data downloaded from arb-silva database
### To the SILVA dataset, we add sequences obtained by assembling ribosomal sequences from Misof et al transcriptomes
### Only the longest sequence is kept per species.
### For each genus, we keep one sequence per nominal species and the longest sequence for unnamed species.
### This script generates files that will be used in further steps:
### Resulting sequences will be placed on Misof et al phylogeny.
### papara is used to extend reference alignments and raxml for phylogenetic placement
### Sequences for which nominal order-level taxonomy does not match their phylogenetic placement are discarded
### This eliminates potential contaminations and phylogenetically uninformative sequences


from Bio import SeqIO, Entrez
import argparse, requests, random, sys, tarfile, glob, signal
sys.path.insert(0,'../')
from get_taxonomy import GNparser, search_name, otl_taxon


# Sometimes NCBIWWW gets stuck for over a day
# Using the signal package we will time it out and raise an exception after 1 hour
# Code adapted from https://stackoverflow.com/questions/25027122/break-the-function-after-certain-time

class TimeOutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeOutError('Call to NCBI timed out!')

signal.signal(signal.SIGALRM, timeout_handler)


#function to get lineage information:
#input can be an ott id, or the response from ott taxon_info
def get_lineage_info(ott_id=None,ott_response=None):
    outdict = dict()
    if ott_response:
        r = ott_response
        ott_id = r.json()['ott_id']
    else:
        r = otl_taxon(ott_id)
        
    if not r:
        return None

    if r.json()['rank'] == 'species':
        outdict['s_ott_id'] = ott_id
    elif r.json()['rank'] == 'genus':
        outdict['genus'] = r.json()['name']
        outdict['g_ott_id'] = ott_id
    elif r.json()['rank'] == 'subgenus' and ('genus' not in [lineage['rank'] for lineage in r.json()['lineage']]): #there is at least one case of a species with a subgenus but no genus. In this case, will consider subgenus as genus
        outdict['genus'] = r.json()['name']
        outdict['g_ott_id'] = ott_id
    elif r.json()['rank'] == 'family':
        outdict['family'] = r.json()['name']
        outdict['f_ott_id'] = ott_id

    for lineage in r.json()['lineage']:
        if lineage['rank'] == 'species':
            outdict['species'] = lineage['name']
            outdict['s_ott_id'] = str(lineage['ott_id'])
        elif lineage['rank'] == 'genus':
            outdict['genus'] = lineage['name']
            outdict['g_ott_id'] = str(lineage['ott_id'])
        elif lineage['rank'] == 'subgenus' and ('genus' not in [lineage['rank'] for lineage in r.json()['lineage']]):#there is at least one case of a species with a subgenus but no genus. In this case, will consider subgenus as genus
            outdict['genus'] = lineage['name']
            outdict['g_ott_id'] = str(lineage['ott_id'])
        elif lineage['rank'] == 'family':
            outdict['family'] = lineage['name']
            outdict['f_ott_id'] = str(lineage['ott_id'])
        elif lineage['rank'] not in ['no rank', 'phylum', 'kingdom', 'domain']:
            outdict[lineage['rank']] = lineage['name']

    return outdict

#this function retrieves data from a fully annotated sequence name to a dictionary
def parse_seqname(name):
    name_split = name.split('|')
    if name_split[0] in ['18S','28S']: #18S and 28S sequences start with a '18S'or '28S'
        name_dict = {'gene':name_split.pop(0),'accession':name_split.pop(0), 'original_name':name_split.pop(-1)}
    elif ':' in name_split[0]: #in reference tree names, there is no accession number
        name_dict = {'original_name':name_split.pop(-1)}
    else: #for most names, first record is simply accession, last is original name, and in the middle there is a series of pairs of key:value
        name_dict = {'accession':name_split.pop(0), 'original_name':name_split.pop(-1)}
    name_dict.update({x.split(':')[0]:x.split(':')[1] for x in name_split})
    return name_dict

if __name__ == "__main__":
    #only mandatory argument is gene (18S or 28S)
    parser = argparse.ArgumentParser()
    parser.add_argument('gene', help = '18S or 28S')
    parser.add_argument('-s','--seed', help = 'random seed', default = 12534)
    parser.add_argument('-g', '--gnparser', help = 'path to GNparser', default = '../gnparser-0.3.1/bin/gnparse')
    parser.add_argument('-e','--email', help = 'email to provide to ncbi')

    #args = parser.parse_args(['18S']) #uncomment for testing
    args = parser.parse_args()
    gnpath = args.gnparser
    Entrez.email = args.email
    
    #open SILVA alignment from SILVA tarfile
    print 'READING SILVA FILE'
    tar = tarfile.open(glob.glob(args.gene + '/*.tgz')[0], 'r:gz')
    #f = tar.extractfile(tar.getmembers()[1]) #uncomment for testing in mac-generated tarfile
    f = tar.extractfile(tar.getmembers()[0])
    silva = [record for record in SeqIO.parse(f,'fasta')]
    f.close()
    tar.close()

    #back transcribe to DNA
    for record in silva:
        record.seq = record.seq.back_transcribe()
    
    print 'FILTERING TO ONE SEQUENCE PER SPECIES'

    #do a first quick filter to one sequence per nominal species to reduce load on open tree API
    accessions = [record.description.split('.')[0] for record in silva]
    names = [record.description.split(';')[-1] for record in silva]
    lengths = [len(record) for record in silva]

    derep_seqs = []
    derep_names = []
    derep_accessions = []
    for name in set(names):
        random.seed(str(args.seed) + name) #to make choice random and yet reproducible, random seed is the seed provided concatenated with organism name
        indexes_to_keep = [i for i in xrange(len(silva)) if name == names[i]]
        maxlen = max([lengths[i] for i in indexes_to_keep])
        indexes_to_keep = [i for i in indexes_to_keep if lengths[i] == maxlen]
        index = random.choice(indexes_to_keep)
        derep_seqs.append(silva[index][:])
        derep_names.append(name)
        derep_accessions.append(accessions[index])
        

    print 'SEARCHING OPEN TREE OF LIFE'
    #save ott version
    ott_version = 'ott_version:' + str(requests.post('https://api.opentreeoflife.org/v3/taxonomy/about').json()['source'])

    #retrieve ncbi id for each accession and then search in ott to get ott id
    seq_info = [dict() for i in xrange(len(derep_seqs))]
    for i in xrange(len(derep_accessions)):
        ott_id = None
        r = None
        ncbi_id = None

        print 'Searching ' + str(derep_accessions[i]) + ' on NCBI'
        signal.alarm(30)
        while True:
            try:
                #first, try getting ott id by the ncbi accession
                try:
                    entrez_id = Entrez.read(Entrez.esearch(db='nucleotide',term= derep_accessions[i]))['IdList']
                    ncbi_id = Entrez.read(Entrez.efetch(db='nucleotide',id=entrez_id, retmode = 'xml', rettype='fasta', seq_stop=1))[0]['TSeq_taxid']
                    r = otl_taxon(ncbi_id, ncbi = True)
                    ott_id = r.json()['ott_id']
                    current_name = r.json()['unique_name']
                    lineage_info = get_lineage_info(ott_response=r)
                    if 's_ott_id' and 'g_ott_id' in lineage_info.keys():
                        seq_info[i].update({'current_name':current_name})
                        seq_info[i].update(lineage_info)
                        print derep_names[i] + ' OK. Found in ott as ' + current_name
                    elif 'g_ott_id' in lineage_info.keys():
                        seq_info[i].update({'current_name':current_name})
                        seq_info[i].update(lineage_info)
                        print derep_names[i] + ' OK. Only genus found in ott as ' + current_name
                    else:
                        print derep_names[i] + ': genus not found in ott, skipping \r'
                        break            
                    
                    
                #if it does not work, use a name search    
                except:
                    print 'accession not found, searching by name: ' + str(derep_names[i])
        
                    try:
                        cannonical_name = GNparser(derep_names[i],gnpath)
                    except KeyError:
                        print derep_names[i] + ': could not parse name, skipping \r'
                        break
            
                    if 'cg' not in cannonical_name.keys():
                        print derep_names[i] + ': could not recognize genus name, skipping \r'
                        break
            
                    try:
                        search_response = search_name(cannonical_name['cg'],cannonical_name['cs'],gnpath,context = 'Arthropods', genus_only = False)
                        if search_response['tax_source'] == 'OTT':
                            matched_name = search_response['matched_name']
                            current_name = search_response['current_name']
                            ott_id = search_response['source_id']
                            provider = search_response['tax_source']
                            taxlevel = search_response['tax_level']
                    except KeyError: #if not 'cs' or species not in OTT, try looking for genus
                        try:
                            search_response = search_name(cannonical_name['cg'],'',gnpath,context = 'Arthropods', genus_only = True)
                            if search_response['tax_source'] == 'OTT':
                                matched_name = search_response['matched_name']
                                current_name = search_response['current_name']
                                ott_id = search_response['source_id']
                                provider = search_response['tax_source']
                                taxlevel = search_response['tax_level']
                        except (TypeError,KeyError): #if response is None, not found; if KeyError, not found either (since 'cg' should exist, error is within search_name function)
                            print derep_names[i] + ': genus could not be found as an arthropod, skipping \r'
                            break
                    except TypeError: #if response is None, not found
                        print derep_names[i] + ': species could not be found as an arthropod, skipping \r'
                        break
            
                    if ott_id: #add lineage info
                        lineage_info = get_lineage_info(ott_id)
                        if taxlevel == 'species':
                            seq_info[i].update({'current_name':current_name})
                            seq_info[i].update(lineage_info)
                            print derep_names[i] + ' OK. Found in ott as ' + current_name
                        elif taxlevel == 'genus':
                            seq_info[i].update({'current_name':current_name})
                            seq_info[i].update(lineage_info)
                            print derep_names[i] + ' OK. Only genus found in ott as ' + current_name
                        else:
                            print derep_names[i] + ': genus not found in ott, skipping \r'
                            break
                    else:
                        print derep_names[i] + ': not found in ott, skipping \r'
                        break
                #remove any keys that are None
                for key in seq_info[i].keys()[:]:
                    if seq_info[i][key] is None:
                        del seq_info[i][key]
                break #this leaves the while loop, in case of no timeout error
            except TimeOutError as err:
                print err
                print 'Trying again'
        signal.alarm(0)

    

    #now delete records not found on ott and rename sequences with taxonomic information
    for i in reversed(xrange(len(derep_seqs))):
        if not seq_info[i]:
            seq_info.pop(i)
            derep_seqs.pop(i)
            derep_names.pop(i)
        else:
            derep_seqs[i].id = '|'.join([derep_seqs[i].name] + sorted([k + ':' + str(x) for k,x in seq_info[i].iteritems()]) + [ott_version,derep_names[i]])
            derep_seqs[i].description = ''
            derep_seqs[i].name = ''

    #load reference sequences if they are note present already and add information about taxonomic ranks
    accessions = [parse_seqname(seq.id)['accession'] for seq in derep_seqs]
    ref_seqs = [record for record in SeqIO.parse(args.gene + '/' + args.gene + '_renamed.fasta','fasta')]
    new_ref_seqs = []
    for record in ref_seqs[:]:
        record.id = record.description
        record.description = ''
        record.name = ''
        record_dict = parse_seqname(record.id)
        if record_dict['accession'] in accessions:
            continue
        record_dict.update(get_lineage_info(record_dict['g_ott_id']))
        gene = record_dict.pop('gene')
        original_name = record_dict.pop('original_name')
        accession = record_dict.pop('accession')
        record.id = '|'.join([accession] + sorted([k + ':' + str(x) for k,x in record_dict.iteritems()]) + [ott_version,original_name])
        new_ref_seqs.append(record)
    ref_seqs = new_ref_seqs
        

    #now we do a second round of name-based filtering, using s_ott_id and g_ott_id
    #this time we also add reference sequences
    #for each g_ott_id, if one or more s_ott_ids, keep one sequence for each (the longest)
    #if no sequences have an s_ott_id, keep the longest one
    all_seqs = derep_seqs + ref_seqs
    final_seqs = []
    genus_dict = dict()
    for record in all_seqs:
        try:
            g_ott,s_ott = [int(parse_seqname(record.id)[key]) for key in ['g_ott_id','s_ott_id']]
        except :
            #print record.id
            g_ott = int(parse_seqname(record.id)['g_ott_id'])
            s_ott = None
        if g_ott in genus_dict.keys():
            genus_dict[g_ott].append([s_ott,record])
        else:
            genus_dict[g_ott] = [[s_ott,record]]

    for g_ott,x in genus_dict.iteritems():
        random.seed(args.seed + g_ott) #to make choice random and yet reproducible, random seed is the seed provided in the script plus genus ott_id
        if all([record[0] is None for record in x]): #if none of species for this genus is identified, keep only the longest sequence
            lengths = [len(record[1]) for record in x]
            index_to_choose = random.choice([i for i in xrange(len(lengths)) if lengths[i] == max(lengths)])
            final_seqs.append(x[index_to_choose][1])
        else: #if some or all of the species for this genus are identified, keep only one sequence por s_ott_id
            s_ott_ids = set([record[0] for record in x if record[0] is not None])
            for s_ott in s_ott_ids:
                temp_list = [record[1] for record in x if record[0] == s_ott]
                lengths = [len(record) for record in temp_list]
                index_to_choose = random.choice([i for i in xrange(len(lengths)) if lengths[i] == max(lengths)])
                final_seqs.append(temp_list[index_to_choose])

    #finally, remove gene name from selected seq names
    for record in final_seqs:
        record_dict = parse_seqname(record.id)
        try:
            del record_dict['gene']
        except KeyError:
            pass
        original_name = record_dict.pop('original_name')
        accession = record_dict.pop('accession')
        record.id = '|'.join([accession] + sorted([k + ':' + str(x) for k,x in record_dict.iteritems()]) + [ott_version,original_name])


    #save sequences to align and to use raxml for an additional filtering based on phylogenetic placement
    SeqIO.write(final_seqs,open(args.gene + '/query_unaligned.fasta','w'),'fasta')

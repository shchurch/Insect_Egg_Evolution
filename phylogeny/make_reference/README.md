# Making a alignments with ribosomal sequences from genera in Misof et al 2015

By Bruno de Medeiros, starting on November 2016


This folder contains scripts to obtain reference sequences for 18S and 28S, and also other genes

## 18S and 28S

This folder has scripts to obtain 18S and 28S ribosomal sequences for species included in the phylogeny of Misof et al. Most of the sequences were obtained by assembling from SRA files indicated in their supplement. Fo species not there, I obtained ribosomal sequences from congeneric species from genbank or by assembling other SRA files.
The following pipeline was used to assemble these sequences.

SLURM scripts were run in Harvard Odyssey Cluster


1. **accessions.csv** has SRA accessions and species names.
2. **download_rRNA.sh** downloads SRA files
3. **assemble_rRNA.sh** uses ncbi program `fastq-dump` to generate fastq files, `Trimmomatic` to clean sequences, `bowtie2` to separate reads matching an rRNA dataset downloaded from SILVA (in folder *bowtie*) and Trinity to assemble reads into contigs.
4. **select_sequence.sh** uses CD-HIT to remove redundant contigs and chooses the longer contig as the sequence for each sample (18S and 28S were assembled separately).
5. resulting sequences are loaded into Geneious, and manually checked by BLASTING. If too short or BLAST indicates contamination, replaced by a sequence from congeneric species obtained from genbank. 
6. final selection of sequences obtained by de novo assembly or searching ncbi is exported from Geneious to folder `final/`. All of the following steps take place within this folder:
7. **download_tree_and_rename_seqs.py** is used to:
  * standardize sequence names, including taxonomic information from open tree taxonomy
  * download tree from misof et al, mapped to open tree taxonomy
  * prune misof tree only to genera present in the reference sequences
  * export pruned trees and sequences with fully annotated sequence and tip names
  * export a version of the trees and sequences with genus ott id as sequence and tip names for papara/raxml
7. Reference 18S and 28S sequences for RAxML placement algorithm are loaded into Geneious, and a a new blast search is conducted to check orientation. Reverse complement is obtained as needed, sequences are re-exported.
  * during the new blast, I discovered that I had overlooked a full-length Stylops 18S from genbank, substituted that with my original assembled sequence
8. **reorient_seqs.py** makes a new set of fully annotated sequences in the correct orientation, based on the sequences exported from Geneious.
9. 28S alignment has a few fragmented sequences, so it was refined by hand. The only manual change was deletion of isolated small blocks of nucleotides at both ends of some sequences.


### List of assembled sequences found to be problematic:
  * 18S (<1200bp):
Acerentomon - blasting to parasitoid wasps
Campodea - blasting to fungus

  * 18S (>1200bp):
Annulipalpia - blasting to Hymenoptera
Atelura - blasting to Odonata
Occasjapyx - blasting to mayflies

  * 28S (<1900 bp):
Acanthocasuarina - blasting to stick insects
Acerentomon - blasting to earwigs
Annulipalpia - blasting to Mecoptera
Anurida - blasting to plants
Campodea - blasting to fungi
Conwentzia - blasting to many insects, uninformative
Corydalus - blasting to many insects, uninformative
Ectopsocus - blasting to many insects, uninformative
Euroleon - blasting to many insects, uninformative
Liposcelis - blasting to many insects, uninformative
Meinertellus - blasting to many insects, uninformative
Menopon - blasting to many insects, uninformative
Occasjapyx - blasting to mayflies
Okanagana - blasting to many insects, uninformative
Osmylus - blasting to many insects, uninformative
Trialeurodes - blasting to moths
Zorotypus - blasting to ants
Bemisia - blasting to plants

  * 28S (>1900 bp):
none

### Sequences excluded from alignment:
  * Annulipalpia - no substitute found on genbank - no problem, many trichoptera
  * Atelura - no substitute found on genbank - no problem, other Zygentoma
  * Acanthocasuarina - no substitute found on genbank - no problem, plenty of hemiptera
  * Trialeurodes - no substitute found on genbank - no problem, plenty of hemiptera
  * Bemisia - no substitute found on genbank
  
### Sequences assembled from SRA, but replaced from genbank and/or SILVA reference sequences:
  * Acerentomon - 18S and 18S
  * Campodea - 18S and 28S
  * Occasjapyx - 18S and 28S
  * Anurida - 28S
  * Zorotypus - 28S
  * Priacma - 18S and 28S
  * Mengenilla - 18S and 28S
  * Carabus - 18S and 28S
  * Stylops - 18S
  * Drosophila - 28S
  * Tricholepidion - 28S
  * Bittacus - 28S

### Sequences maybe uninformative that were kept for a first pass:
  * Conwentzia - no substitute found on genbank 
  * Corydalus - no substitute found on genbank 
  * Ectopsocus - no substitute found on genbank 
  * Euroleon - no substitute found on genbank 
  * Liposcelis - no substitute found on genbank 
  * Meinertellus - no substitute found on genbank 
  * Menopon - no substitute found on genbank 
  * Okanagana - no substitute found on genbank 
  * Osmylus - no substitute found on genbank

# Genus-level phylogeny

By Bruno de Medeiros, starting on November 2016

This folder contains scripts used to generate a phylogeny for insect genera present in the insect egg database.


## Genes

We used 18S and 28S rRNA sequences downloaded from the [Silva ribosomal database](https://www.arb-silva.de), a curated database of rRNA sequences. For both ribosomal subunits, we browsed sequences on Silva website and downloaded sequences using the following filters:
* RefNR
* taxonomy: **Hexapoda**
* length: **>500**
* alignment quality: **>50**


## Folder contents

This folder contains the following sub folders:
* **make_reference:**  scripts and files to produce a reference dataset to classify sequences downloaded from Silva. Check README file in the folder for details.
* **raxml_epa:** scripts used for phylogenetic placement using the references produced in *make_reference*. 
* **18S** contains downloaded sequences from SILVA and final 18S alignment.
* **28S** same as above, for 28S.
* **constraint_info** contains files needed for setting up phylogenetic constraints based on the Misof and Rainford trees
* **mrbayes_smalltrees** contains scripts to run mrbayes
* **phylo_figure** scripts to make phylogenetic trees in the supplement

## Steps

1. **1_silva_to_alignment.py** runs **silva_to_alignment.py** to filter sequences by taxonomy: only one sequence per species ott_id, only one per genus ott_id if no species ott_id.

2. **2_raxml_epa.sh** makes an alignment using MAFFT, adding new sequences to the reference alignment, and then uses RAxML for phylogenetic placement.

3. **3_final_alignments.sh** takes RAxML results and filters out sequences for which phylogenetic placement does not match OTT taxonomy at order level. Alignments are re-built in this smaller dataset using mafft, without Gblocks.
  * this script makes use of mafft and of a python script **get_final_dataset.py**

4. **4_trim_alignments.sh** calls **trim_alignments.py** and trims the ends of alignments to the area with at least 20% species. Within the area kept, it removes sites with less than 10 sequences with data. Finally, it removes sequences with less than 100 bp of non-missing data after trimming.

5. **5_build_raxml_and_mrbayes.sh** builds input files for RAxML and mrbayes. We ended up using only Mrbayes. It first runs **get_concat_alignment_for_phenotypes.py**, which reads final alignments and compares to the egg dataset. Overlapping species and genera, and parittion files, are written to  folder *ML_tree*. This script also writes csv files with species overlap. Following, it runs **build_mrbayes.py**, which builds files to run mrbayes using genus-level data. It reads alignments, plus Rainford family-level tree and Misof order-level tree and makes a mrbayes input file including constraints from Rainford (hard for families monophyletic in OTL, partial between monophyletic families, hard for orders and above). Also calibrates the tree with between-order estimates from rainford or Misof.

7. **mrbayes_smalltrees/get_mrbayes_trees_4runs.sh** bash script that runs mrbayes on Odyssey cluster and automatically resumes if run stops due to timeout. Before resuming, it checks parameter and tree files for errors and deletes the last samples of MCMC if needed.

8. **6_stitch_mrbayes.sh** gets small trees from folder **mrbayes_smalltrees** and grafts them onto the appropriate backbone by using **stitch_mrbayes.py**

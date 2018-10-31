This folder contains files used for constraints and calibrations

* misof_full_annotated.tre

Misof maximum likelihood tree downloaded from open tree of life with tips renamed to retain all taonomic information. Tip names can be transformed to a dictionary by using the function parse_seqname() in silva_to_alignment.py

See folder make_reference for scripts to obtain this tree.

* misof_figure_tree_topology.tsv

File with the tree topology manually obtained from Misof et al figure, used to put node numbers and time estimates in misof_full_annotated.tre


* misof_name_translation.csv

File to help relating names in misof_full_annotated.tre to those in misof et al's figure

* misof_File_S17.txt

Supplementary file S17 from Misof et all, with estimated node ages

* misof_full_with_nodes.nexml

Tree obtained by combining the information of all the above.

Nodes are named, and estimated ages are added as annotations

Tips are named as genus_ottid, and all additional information is added as annotations


* annotated_misof_tree.py 

Script that produces misof_full_with_nodes.nexml


* rainford_family_tree.nxs

Calibrated family-level tree from Rainford et al

* rainford_family_key.txt

Translation of Rainford's tip names to full family/order names.


# NAR_splicing

### Scripts used in "Structural disruption of exonic stem-loops immediately upstream of the intron regulates mammalian splicing"

sjshapePC.py collects icSHAPE enrichment scores upstream of splice junctions.
intronHunterPC.py collects icSHAPE enrichment scores upstream of annotated retained introns.

the *_3prime.py versions do the same, but for downstream regions.

Usage is the same for all 4 scripts: `scriptname.py icSHAPE_file gtf_file number_of_threads output_file`.

where `icSHAPE_file` is the final output from the [icSHAPE pipeline](https://github.com/qczhang/icSHAPE).

flankPU.py takes the output of one of the above scripts, extracts the sequence plus up to a specified length of flanking sequence up- and downstream, and uses [RNAplfold](https://github.com/ViennaRNA/ViennaRNA) to determine PU values for all hexamers across the region of interest.
Works for 5' or 3' regions; last argument must specify "5" or "3".

Usage is `flankPU.py region_file fasta_file gtf_file output_file_base flank_length 5|3`

where `region_file` is an output file from sjshapePC* or intronHunterPC*.

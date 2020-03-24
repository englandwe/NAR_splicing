# NAR_splicing

### Scripts used in "Structural disruption of exonic stem-loops immediately upstream of the intron regulates mammalian splicing"

sjshapePC.py collects icSHAPE enrichment scores upstream of splice junctions.
intronHunterPC.py collects icSHAPE enrichment scores upstream of annotated retained introns.

the *_3prime.py versions do the same, but for downstream regions.

Usage is the same for all 4 scripts: `scriptname.py icSHAPE_file gtf_file number_of_threads output_file`.

where `icSHAPE_file` is the final output from the [icSHAPE pipeline](https://github.com/qczhang/icSHAPE).

The flankFetch_2way* scripts take the output of one of the above scripts, extract the sequence plus up to a specified length of flanking sequence up- and downstream, and use one of three methods to determine PU values for all hexamers across the region of interest.

flankFetch_2way_rnaplfold.py uses [RNAplfold](https://github.com/ViennaRNA/ViennaRNA).
flankFetch_2way_localfold.py uses [LocalFold](https://github.com/BackofenLab/LocalFold).

Works for 5' or 3' regions; last argument must specify "5" or "3".

Usage is `flankFetch_2way_*.py region_file fasta_file gtf_file output_file_base flank_length 5|3`

where `region_file` is an output file from sjshapePC* or intronHunterPC*.

flankFetch_2way_memeris.py uses a modified GetSecondaryStructureValues from [MEMERIS](https://github.com/BackofenLab/MEMERIS) to calculate PU values for the desired sequence. This is part of a three-step process:

1. for each desired flank length: flankFetch_2way_memeris.py (usage as above)
2. for each fasta file produced by step 1: perl -w ../GetSecondaryStructureValues_fix.perl -f file.fa -o file.pu -l 6 -method PU
3. meanPU_memeris.py

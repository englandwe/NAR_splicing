# NAR_splicing

### Scripts used in "Structural disruption of exonic stem-loops immediately upstream of the intron regulates mammalian splicing"

sjshapePC.py collects icSHAPE enrichment scores upstream of splice junctions.
intronHunterPC.py collects icSHAPE enrichment scores upstream of annotated retained introns.

the *_3prime.py versions do the same, but for downstream regions.

Usage is the same for all scripts: `scriptname.py icSHAPE_file gtf_file number_of_threads output_file`.

where `icSHAPE_file` is the final output from the [icSHAPE pipeline](https://github.com/qczhang/icSHAPE).

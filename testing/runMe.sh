#!/bin/bash

set -ex

# define functional regions
../define_functional_regions_of_genome.py  --input_spec input_spec.tsv

# build minigenome fasta
../build_minigenome_from_intervals.py  --intervals merged_intervals.tsv --genome genome.fa

# translate original ref annotation coords to minigenome:
../translate_fullgenome_to_minigenome_annot.py --fullgenome_annot ref_annot.gtf --translation_intervals minigenome.coord_translation.tsv --output_gtf minigenome.gtf 


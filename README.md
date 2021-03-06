# ShrinkGenome
Shrink a large genome down to predefined functional regions


Involves three steps:

1.  Defining the functional regions:

```
./define_functional_regions_of_genome.py -h

usage: define_functional_regions_of_genome.py [-h] --input_spec INPUT_SPEC
                                              [--output OUTPUT]

define functional regions of the genome

optional arguments:
  -h, --help            show this help message and exit
  --input_spec INPUT_SPEC
                        file containing input specification with format:
                        interval_filename flank_extend min_feature_length
                        (default: None)
  --output OUTPUT, -O OUTPUT
                        output intervals file (default: merged_intervals.tsv)
```

Interval files contain simple interval files with fields (and 1-based coordinates):
Chromosome   Start   End



>The intervals that define functional reagions are user-defined.  Typically they would involve reference annotations (with settings for 1kb flank extension, and zero minimum feature lengths), coupled to other intervals derived from expression-defined genomic regions or other.  Just don't leave out any region that will be useful to you! :-)


The above generates an output file containing the merged intervals (eg. merged_intervals.tsv)


2.  Build the mini-genome according to the functional regions as defined.

```
./build_minigenome_from_intervals.py -h

usage: build_minigenome_from_intervals.py [-h] --intervals INTERVALS --genome
                                          GENOME [--spacer_len SPACER_LEN]

build minigenome fasta

optional arguments:
  -h, --help            show this help message and exit
  --intervals INTERVALS
                        intervals bed file (default: None)
  --genome GENOME       genome fasta file (default: None)
  --spacer_len SPACER_LEN
                        length of N spacer between sequence regions (default:
                        10)

```


The above generates files:
-  minigenome.fa  : the mini genome fasta file
-  minigenome.coord_translation.tsv  : provides coordinate translations between the original genome and the minigenome.


3. Translate reference genome annotations to the minigenome

Use the following to convert your reference gtf file to the minigenome version.  Note, ideally, this ref genome file is restricted to coordinates that were included when defining the functional regions for which the minigenome was defined.

```
./translate_fullgenome_to_minigenome_annot.py -h

usage: translate_fullgenome_to_minigenome_annot.py [-h] --fullgenome_annot
                                                   FULLGENOME_ANNOT
                                                   --translation_intervals
                                                   TRANSLATION_INTERVALS
                                                   --output_gtf OUTPUT_GTF

translate full genome gtf|gff3 to the minigenome gtf|gff3

optional arguments:
  -h, --help            show this help message and exit
  --fullgenome_annot FULLGENOME_ANNOT
                        input fullgenome gtf or gff3 file to translate coords
                        (default: None)
  --translation_intervals TRANSLATION_INTERVALS
                        translation intervals tsv file (default: None)
  --output_gtf OUTPUT_GTF
                        output gtf filename (default: None)

```



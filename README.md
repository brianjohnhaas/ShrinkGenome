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
usage: build_minigenome_from_intervals.py [-h] --intervals_bed INTERVALS_BED
                                          --genome GENOME
                                          [--spacer_len SPACER_LEN]

build minigenome fasta

optional arguments:
  -h, --help            show this help message and exit
  --intervals_bed INTERVALS_BED
                        intervals bed file (default: None)
  --genome GENOME       genome fasta file (default: None)
  --spacer_len SPACER_LEN
                        length of N spacer between sequence regions (default:
                        10)

```




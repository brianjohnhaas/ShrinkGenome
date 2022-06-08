#!/usr/bin/env python3

import sys, os, re
import subprocess
import logging
import argparse
import pandas as pd
import pyranges as pr

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="build minigenome fasta",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--intervals_bed", type=str, required=True, help="intervals bed file")
    parser.add_argument("--genome", type=str, required=True, help="genome fasta file")
    parser.add_argument("--spacer_len", type=int, required=False, default=10, help="length of N spacer between sequence regions")

    args = parser.parse_args()

    intervals_bed_file = args.intervals_bed
    genome_fasta_file = args.genome
    Nspacer_length = args.spacer_len

    if not os.path.exists(f"{genome_fasta_file}.fai"):
        logger.info(f"Creating fasta index {genome_fasta_file}.fai")
        subprocess.check_call(f"samtools faidx {genome_fasta_file}", shell=True)


    logger.info("parsing {}".format(intervals_bed_file))
    intervals = pd.read_csv(intervals_bed_file, sep="\t", header=None)
    intervals = intervals.iloc[:,[0,1,2]]
    intervals.columns = ["Chromosome", "Start", "End"]

    intervals = intervals.astype({'Start':'int', 'End':'int'})
    
    # just be sure they're non-overlapping
    logger.info("merging intervals just in case...")
    pr_intervals = pr.PyRanges(intervals).merge(strand=False)
    intervals = pr_intervals.df.copy()
    pr_intervals = None

    logger.info("sorting intervals")
    intervals = intervals.sort_values(['Chromosome', 'Start'])

    chromosomes = intervals.Chromosome.unique()
    logger.info("Have chromosomes: {}".format(chromosomes))
    
    Nspacer = "N" * Nspacer_length

    minigenome_bed_ofh = open("minigenome.coord_translation.tsv", 'wt')
    print("\t".join(["Chromosome", "Start", "End", "NewStart", "NewEnd"]), file=minigenome_bed_ofh)
    
    minigenome_fa_ofh = open("minigenome.fa", "wt")
    
    for chromosome in chromosomes:

        logger.info("-processing {}".format(chromosome))
        
        chr_data = intervals[intervals.Chromosome == chromosome]

        contig_new_lend = 0
        contig_new_rend = 0

        print(f">{chromosome}", file=minigenome_fa_ofh)

        chromosome_seq = extract_genome_seq(genome_fasta_file, chromosome)
        
        for i, row in chr_data.iterrows():

            lend = max(row.Start, 1) # avoid any incoming negative coordinates
            rend = row.End
            segment_len = rend - lend + 1
            
            contig_new_lend = contig_new_rend + 1
            contig_new_rend = contig_new_lend + segment_len - 1
            
            print("\t".join([chromosome, str(lend), str(rend), str(contig_new_lend), str(contig_new_rend)]), file=minigenome_bed_ofh)

            # add short N spacer
            
            genome_seq = chromosome_seq[lend-1:rend]  # remember indexing of sequence is from zero!
            genome_seq += Nspacer

            contig_new_rend += Nspacer_length
            
            print(genome_seq, file=minigenome_fa_ofh, end='')

        print("", file=minigenome_fa_ofh) # newline seq separator


    logger.info("done")
        
    sys.exit(0)



def extract_genome_seq(genome_fasta_file, chromosome):

    cmd = f"samtools faidx {genome_fasta_file} {chromosome}"
    logger.info(cmd)
    fa_record = subprocess.check_output(cmd, shell=True).decode()
    fa_record = fa_record.split("\n")
    fa_record = fa_record[1:]
    fa_record = "".join(fa_record)

    return fa_record




if __name__=='__main__':
    main()

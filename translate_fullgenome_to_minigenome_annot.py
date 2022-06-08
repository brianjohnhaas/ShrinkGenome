#!/usr/bin/env python3

import sys, os, re
import subprocess
import logging
import argparse
import pandas as pd
import pyranges as pr
import csv

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="translate full genome gtf|gff3 to the minigenome gtf|gff3",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument("--fullgenome_annot", type=str, required=True, help="input fullgenome gtf or gff3 file to translate coords")
    parser.add_argument("--translation_intervals", type=str, required=True, help="translation intervals tsv file")
    parser.add_argument("--output_gtf", type=str, required=True, help="output gtf filename")
    
    args = parser.parse_args()

    fullgenome_annot_filename = args.fullgenome_annot
    translation_intervals_filename = args.translation_intervals
    output_gtf = args.output_gtf
    
    logger.info("parsing {}".format(fullgenome_annot_filename))

    df = pd.read_csv(fullgenome_annot_filename, sep="\t", names=["Chromosome", "Source", "Type", "Start", "End", "Something", "Strand", "Dot", "Info"])

    logger.info("parsing {}".format(translation_intervals_filename))
    pd_translation = pd.read_csv(translation_intervals_filename, sep="\t") # already has the Chromosome, Start, End headings


    # restrict to those annotation features that are on coordinate-transfer contigs.
    chromosomes_can_be_translated = pd_translation['Chromosome'].unique().tolist()

    pd_translation.set_index('Chromosome', drop=False, inplace=True)

    df = df[ df['Chromosome'].isin(chromosomes_can_be_translated) ]
    

    def overlap_by_contig(group_df):
        #print(group_df)
        #print(group_df.shape)
        chromosome = group_df['Chromosome'].unique().tolist()[0]

        print(chromosome)
        
        translation_df = pd_translation.loc[[chromosome]]
        #print(translation_df)
        
        pr_group_df = pr.PyRanges(group_df)
        pr_translation = pr.PyRanges(translation_df)

        pr_join = pr_group_df.join(pr_translation)

        return pr_join.df.copy()
    
    
    logger.info("joining based on feature coordinate overlaps by chromosome")
    data = df.groupby('Chromosome', as_index=False).apply(overlap_by_contig)

    logger.info("assigning adjusted coordinates")
    data['Start2'] = data['Start'] - data['Start_b'] + data['NewStart']
    data['End2'] = data['End'] - data['Start_b'] + data['NewStart']


    logger.info("adding original coordinates as an annotation")
    
    data['OrigStart'] = data['Start'].astype(int)
    data['OrigEnd'] = data['End'].astype(int)
    
    def add_orig_coords(row):
        return row['Info'] + f" OrigCoords \"{row['OrigStart']}-{row['OrigEnd']}\" "
    
    data['Info'] = data.apply(add_orig_coords, axis=1)
    
        
    
    data['Start'] = data['Start2'].astype(int)
    data['End'] = data['End2'].astype(int)
    data = data[ ["Chromosome", "Source", "Type", "Start", "End", "Something", "Strand", "Dot", "Info"] ]

    logger.info("writing output gtf: {}".format(output_gtf))
    data.to_csv(output_gtf, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
    

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

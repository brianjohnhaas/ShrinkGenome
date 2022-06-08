#!/usr/bin/env python3


import sys, os, re
import pyranges as pr
import pandas as pd
import logging
import argparse

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)




def main():

    parser = argparse.ArgumentParser(description="define functional regions of the genome",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--ref_annot", type=str, required=False, help="reference annotation bed file")
    parser.add_argument("--trans_bulk", type=str, required=False, help="bulk RNA-seq transcripts")
    parser.add_argument("--trans_SC", type=str, required=False, help="single cell RNA-seq transcripts")
    parser.add_argument("--output_bed", "-O", type=str, default="merged.bed", help="output bed file")
    
    args = parser.parse_args()

    ref_annot_bed = args.ref_annot
    trans_bulk_bed = args.trans_bulk
    trans_SC_bed = args.trans_SC
    output_bed_filename = args.output_bed

    if not (ref_annot_bed or trans_bulk_bed or trans_SC_bed):
        exit("use --help for menu")
    

    ref_annot_extend = 1000

    trans_bulk_min_interval_size = 100
    trans_bulk_extend = 100

    trans_SC_min_interval_size = 100
    trans_SC_extend = 100

    pr_merged = None

    if ref_annot_bed:
        pr_merged = parse_merged_intervals(ref_annot_bed, ref_annot_extend, 0)

    if trans_bulk_bed:
        pr_trans_bulk = parse_merged_intervals(trans_bulk_bed, trans_bulk_extend, trans_bulk_min_interval_size)
        if pr_merged:
            pr_merged = join_pair_of_intervals(pr_merged, pr_trans_bulk)
        else:
            pr_merged = pr_trans_bulk


    if trans_SC_bed:
        pr_trans_SC = parse_merged_intervals(trans_SC_bed, trans_SC_extend, trans_SC_min_interval_size)
        if pr_merged:
            pr_merged = join_pair_of_intervals(pr_merged, pr_trans_SC)
        else:
            pr_merged = pr_trans_SC

    logger.info("-Writing output intervals")
    pr_merged.to_bed(output_bed_filename)
    logger.info("-done")

    sys.exit(0)


def parse_merged_intervals(bed_file, extend=0, min_interval_size=0):

    logger.info("Parsing {}".format(bed_file))
    data = pd.read_csv(bed_file, sep="\t", names=["Chromosome", "Start", "End"])

    if min_interval_size > 0:
        logger.info("-pruning intervals < {} in length".format(min_interval_size))
        data = data[data.End - data.Start >= min_interval_size]

    pr_merged = pr.PyRanges(data)
    pr_merged = pr_merged.merge(strand=False)
    logger.info("Parsed {} with  {:,} merged intervals totaling {:,} bases".format(bed_file, len(pr_merged), pr_merged.length))

    if extend > 0:
        logger.info("Extending intervals by {} on each side".format(extend))
        pr_merged.Start -= extend
        pr_merged.End += extend
        pr_merged = pr_merged.merge(strand=False)
        logger.info("    Now {:,} merged intervals totaling {:,} bases".format(len(pr_merged), pr_merged.length))
    
    return pr_merged


def join_pair_of_intervals(pr_1, pr_2):
    logger.info("-Merging with previous interval set")
    data = pd.concat([pr_1.df, pr_2.df])
    pr_both = pr.PyRanges(data)
    pr_both = pr_both.merge(strand=False)
    logger.info("    Merged intervals count: {:,}  totaling {:,} bases".format(len(pr_both), pr_both.length))
    
    return pr_both


if __name__=='__main__':
        main()

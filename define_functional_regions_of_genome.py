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

    parser.add_argument("--input_spec", type=str, required=True,
                        help="file containing input specification with format:\n" +
                        "interval_filename\tflank_extend\tmin_feature_length\n\n")

    parser.add_argument("--output", "-O", type=str, default="merged_intervals.tsv", help="output intervals file")
    
    args = parser.parse_args()

    input_spec_filename = args.input_spec
    output_filename = args.output
    

    targets = list()
    with open(input_spec_filename) as fh:
        for line in fh:
            line = line.rstrip()
            
            if re.match("#", line):
                continue
            
            if re.match(r"\w", line) is None:
                continue
            
            filename, flank_extend, min_feature_length = line.split("\t")
            targets.append([filename, int(flank_extend), int(min_feature_length)])


    # define functional regions:        

    pr_merged = None

    for target in targets:
        filename, flank_extend, min_feature_length = target
        
        pr_entry = parse_merged_intervals(filename, flank_extend, min_feature_length)
        
        
        if pr_merged:
            pr_merged = join_pair_of_intervals(pr_merged, pr_entry)
        else:
            pr_merged = pr_entry


    logger.info("-Writing output intervals")
    pr_merged.to_csv(output_filename, sep="\t", header=False)
    logger.info("-done")

    sys.exit(0)


def parse_merged_intervals(intervals_file, extend=0, min_interval_size=0):

    logger.info("Parsing {}".format(intervals_file))
    data = pd.read_csv(intervals_file, sep="\t", names=["Chromosome", "Start", "End"])

    if min_interval_size > 0:
        logger.info("-pruning intervals < {} in length".format(min_interval_size))
        data = data[data.End - data.Start >= min_interval_size]

    pr_merged = pr.PyRanges(data)
    pr_merged = pr_merged.merge(strand=False)
    logger.info("Parsed {} with  {:,} merged intervals totaling {:,} bases".format(intervals_file, len(pr_merged), pr_merged.length))

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

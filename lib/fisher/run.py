#! /usr/bin/env python

import sys
import os
import argparse
import logging
import fisher

def run_compare(args):

    level = logging.getLevelName( args.log_level )

    if args.log_file:
        logging.basicConfig( filename   = args.log_file,
                             level      = level,
                             format     = '%(asctime)s %(message)s',
                             datefmt    ='%m/%d/%Y %I:%M:%S%p' )
    else:
        logging.basicConfig( level      = level,
                             format     = '%(asctime)s %(message)s',
                             datefmt    ='%m/%d/%Y %I:%M:%S%p' )
    #
    # Main function
    #
    fisher.Pileup_and_count( 
            in_bam1 = args.bam1,
            in_bam2 = args.bam2,
            out_file = args.output,
            ref_fa = args.ref_fa,
            baseq_thres = args.base_quality,
            mismatch_rate_disease = args.min_allele_freq,
            mismatch_rate_normal = args.max_allele_freq,
            post_10_q = None,
            fisher_threshold = args.fisher_value,
            min_depth = args.min_depth,
            print_header = args.print_header,
            mapq_thres = args.mapping_quality,
            min_variant_read = args.min_variant_read,
            samtools = args.samtools_path,
            region = args.region
          )


def run_single(args):

    level = logging.getLevelName( args.log_level )

    if args.log_file:
        logging.basicConfig( filename   = args.log_file,
                             level      = level,
                             format     = '%(asctime)s %(message)s',
                             datefmt    ='%m/%d/%Y %I:%M:%S%p' )
    else:
        logging.basicConfig( level      = level,
                             format     = '%(asctime)s %(message)s',
                             datefmt    ='%m/%d/%Y %I:%M:%S%p' )
    #
    # Main function
    #
    fisher.Pileup_and_count( 
            in_bam1 = args.bam1,
            in_bam2 = None,
            out_file = args.output,
            ref_fa = args.ref_fa,
            baseq_thres = args.base_quality,
            mismatch_rate_disease = args.min_allele_freq,
            mismatch_rate_normal = None,
            post_10_q = args.post_10_q,
            fisher_threshold = None,
            min_depth = args.min_depth,
            print_header = args.print_header,
            mapq_thres = args.mapping_quality,
            min_variant_read = args.min_variant_read,
            samtools = args.samtools_path,
            region = args.region
          )


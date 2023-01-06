#! /usr/bin/env python

import sys
import os
import argparse
import logging
from . import fisher

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

    is_anno = True if args.print_format == 'anno' else False
    if not is_anno:
        if args.sample1 == None or args.sample2 == None:
            raise ValueError('--sample1 and --sample2 are required for vcf.') 

    if args.regions != None and args.positions != None:
        raise ValueError('--regions and --positions cannot be specified at the same time.') 

    #
    # Main function
    #
    fisher.Pileup_and_count(
            in_bam1 = args.bam1,
            in_bam2 = args.bam2,
            sample1 = args.sample1,
            sample2 = args.sample2,
            out_file = args.output,
            ref_fa = args.ref_fa,
            baseq_thres = args.base_quality,
            mismatch_rate_disease = args.min_allele_freq,
            mismatch_rate_normal = args.max_allele_freq,
            post_10_q = None,
            fisher_threshold = args.fisher_value,
            min_depth = args.min_depth,
            header_flag = args.print_header,
            min_variant_read = args.min_variant_read,
            samtools = args.samtools_path,
            samtools_params = args.samtools_params,
            region = args.region,
            region_file = args.regions,
            positions_bed = args.positions,
            is_anno = is_anno,
            flag_mis_base_0 = args.flag_mis_base_0,
            mismatch_rate_base_0 = args.mismatch_rate_base_0,
            fisher_threshold_log_base_0 = args.fisher_threshold_log_base_0
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

    is_anno = True if args.print_format == 'anno' else False
    if not is_anno:
        if args.sample1 == None:
            raise ValueError('--sample1 is required for vcf.') 

    if args.regions != None and args.positions != None:
        raise ValueError('--regions and --positions cannot be specified at the same time.') 

    #
    # Main function
    #
    fisher.Pileup_and_count(
            in_bam1 = args.bam1,
            in_bam2 = None,
            sample1 = args.sample1,
            sample2 = None,
            out_file = args.output,
            ref_fa = args.ref_fa,
            baseq_thres = args.base_quality,
            mismatch_rate_disease = args.min_allele_freq,
            mismatch_rate_normal = None,
            post_10_q = args.post_10_q,
            fisher_threshold = None,
            min_depth = args.min_depth,
            header_flag = args.print_header,
            min_variant_read = args.min_variant_read,
            samtools = args.samtools_path,
            samtools_params = args.samtools_params,
            region = args.region,
            region_file = args.regions,
            positions_bed = args.positions,
            is_anno = is_anno,
            flag_mis_base_0 = args.flag_mis_base_0,
            mismatch_rate_base_0 = args.mismatch_rate_base_0,
            fisher_threshold_log_base_0 = args.fisher_threshold_log_base_0
          )

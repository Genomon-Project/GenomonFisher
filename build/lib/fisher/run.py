#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
import fisher

############################################################
#
# Main
#
def main(args):

    #
    # logging setup
    #
    # Level     function            value    description
    # CRITICAL  logging.critical()  50      Output only critical errors
    # ERROR     logging.error()     40      Output errors
    # WARNING   logging.warning()   30      Output warnings
    # INFO      logging.info()      20      Output information
    # DEBUG     logging.debug()     10      Output debug information
    # NOTSET                        0       Output all
    #
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
            input_mpileup = args.input_mpileup,
            ref_fa = args.ref_fa,
            threshold = args.base_quality,
            mismatch_rate = args.mismatch_rate,
            post_10_q = args.post_10_q,
            fisher_threshold = args.fisher_value,
            min_depth = args.min_depth,
            compare = args.compare,
            print_header = args.print_header
          )


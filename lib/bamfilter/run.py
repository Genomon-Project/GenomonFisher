#! /usr/bin/env python

import sys
import os
import math
import argparse
import logging
import bamfilter

############################################################
#
# Main
#
def main(args):  

    try:
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
        # Read bam file
        #
        logging.info( "Read bam/sam file." )

        if ( not os.path.exists( args.bam ) ):
            logging.error( "Input file: {file} not found.".format( file = args.bam ) )
            raise

        bamfilter.FilterBam( args.bam, args.out_bam, args.max_indel, args.max_distance, args.map_quality )

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        raise



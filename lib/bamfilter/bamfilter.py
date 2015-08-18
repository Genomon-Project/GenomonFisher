#!/usr/local/package/python2.7/2.7.2/bin/python
"""

This is a template of Python main module.


"""

import sys
import os
import re
from datetime import datetime
import argparse
import logging
import pysam

#
# Global
#
global myself
myself = None


# CIGAR ID
BAM_CMATCH      = 0 # M
BAM_CINS        = 1 # I
BAM_CDEL        = 2 # D
BAM_CREF_SKIP   = 3 # N
BAM_CSOFT_CLIP  = 4 # S
BAM_CHARD_CLIP  = 5 # H
BAM_CPAD        = 6 # P
BAM_CEQUAL      = 7 # =
BAM_CDIFF       = 8 # X

#
# class definitions
#
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
        return value


############################################################
def FilterBam(
        in_bam,
        out_bam,
        max_indel,
        max_distance,
        map_quality
        ):
    """
    Filter bam and save in bam file

    """

    try:
        if not os.path.exists( "{file}.bai".format( file = in_bam ) ):
            pysam.index( in_bam )

        samfile = pysam.Samfile( in_bam, "rb" )

        #
        # 1. Get header
        #
        header = samfile.header
        outbam = pysam.Samfile( out_bam, "wb", header=header )

        #
        # 2. Read sam/bam file and create filtered reads
        #

        # A. Get edit distance to the reference
        # B. Get the number of deletions and insertions.
        # C. Check the sum of deletions and insertions is less than max_indel ( default: 2 ).
        # D. If C. is met and ( the edit distance - bases of deletion and insertions < max_distance )
        #   save the read.
        #
        logging.info( "Fetch reads and filter" )
        for read in samfile.fetch():

            if read.mapq < map_quality:
                next

            tagNM = 0
            if read.tags:
                for tag, value in read.tags:
                    if tag == 'NM':
                        tagNM = value

            size_D = 0
            size_I = 0
            num_del = 0
            num_ins = 0
            if read.cigar:
                for cigar, value in read.cigar:
                    if BAM_CDEL == cigar:
                        num_del += 1
                        if ( not size_D or size_D < value ):
                            size_D = value

                    elif BAM_CINS == cigar :
                        num_ins += 1
                        if ( not size_I or size_I < value ):
                            size_I = value

                if num_del + num_ins < max_indel:
                    if tagNM - size_D - size_I < max_distance:
                        outbam.write( read )

        outbam.close()
        samfile.close()

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        raise



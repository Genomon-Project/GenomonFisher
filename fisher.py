#!/usr/bin/python
"""

bamfilter.py


"""
import sys
import os
import re
import pysam
import scipy.special
from scipy.stats import fisher_exact as fisher
import argparse
import logging
import math

#
# Globals
#
arg = None
target = None
remove_chr = None
filter_quals = None

#
# Class definitions
#

############################################################
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)

        except KeyError:
            value = self[item] = type(self)()

        return value

#
# Subroutines
#

# data_pair IDs
POS_CHR = 0
POS_COORD = 1
POS_REF = 2
POS_DATA1 = 3
POS_DATA2 = 4
POS_FISHER_SNV = 5
POS_FISHER_INS = 6
POS_FISHER_DEL = 7
POS_COUNT = 8

############################################################
def Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth, compare ):

    #
    # mpileup format
    #
    # chr1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
    #
    # 0 chromosome,
    # 1 1-based coordinate,
    # 2 reference base,
    # 3 the number of reads covering the site (1)
    # 4 read bases (1)
    # 5 base qualities (1)
    # 6 the number of reads covering the site (2)
    # 7 read bases (2)
    # 8 base qualities (2)
    #
    global target
    global remove_chr
    global filter_quals

    try:
        #
        # Prepare mpileup data
        #
        mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
        mp_list_len = len( mp_list )
        ref_base_U = mp_list[ 2 ].upper()
        coordinate = mp_list[ 0:3 ]

        #
        # skip if depth is 0
        #
        if mp_list[ 3 ] == '0' or ( mp_list_len > 6 and mp_list[ 6 ] == '0' ):
            return None

        #
        # data_pair IDs
        # POS_CHR = 0
        # POS_COORD = 1
        # POS_REF = 2
        # POS_DATA1 = 3
        # POS_DATA2 = 4
        # POS_FISHER_SNV = 5
        # POS_FISHER_INS = 6
        # POS_FISHER_DEL = 7
        # POS_COUNT = 8
        #
        data_pair = [ mp_list[ 0 ],
                      int( mp_list[ 1 ] ),
                      mp_list[ 2 ],
                      { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'indel': AutoVivification() },
                      { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'indel': AutoVivification() },
                      1.0,
                      'N:1.0',
                      'N:1.0',
                      0 ]


        #
        # Loop for 2 bam file case
        #
        if compare:
            data_pair[ POS_COUNT ] = 2
            input_list = [ ( POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ),
                           ( POS_DATA2, mp_list[ 6 ], mp_list[ 7 ], mp_list[ 8 ] ) ]
        else:
            data_pair[ POS_COUNT ] = 1
            input_list = [ ( POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ) ]

        #
        # position id,
        # mpileup output 4th row(number of read covering the site),
        # 5th row(read bases),
        # 6th row(base quality)
        #
        for data_id, depth, read_bases, qual_list in input_list:

            indel = AutoVivification()

            #
            # Look for deletion/insertion and save info in 'indel' dictionary
            #
            #   ([\+\-])[0-9]+[ACGTNacgtn]+
            #
            # m.group(1): + or - (deletion/insertion)
            # m.group(2): number of deletion/insertion
            # m.group(3): nucleotides
            #
            deleted = 0
            iter = target.finditer( read_bases )
            for m in iter:
                site = m.start()
                type = m.group( 1 )
                num = m.group( 2 )
                bases = m.group( 3 )[ 0:int( num ) ]
                if bases.islower():
                    strand = ( '-', '+' )
                else:
                    strand = ( '+', '-' )

                key = '\t'.join( coordinate + [ bases.upper() ] )
                if type in indel and key in indel[ type ]:
                    indel[ type ][ key ][ strand[ 0 ] ] += 1
                else:
                    indel[ type ][ key ][ strand[ 0 ] ] = 1
                    indel[ type ][ key ][ strand[ 1 ] ] = 0

                read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int(num) + len( num ) + 1 - deleted: ]
                deleted += 1 + len( num ) + int( num )

            #
            # Remove '^.' and '$'
            #
            read_bases = remove_chr.sub( '', read_bases )
            read_bases = read_bases.translate( None, '$' ) 

            #
            # Error check
            #
            if len( read_bases ) != len( qual_list ):
                logging.error( "mpileup data is not good: {0}, {1}".format( mpileup, read_bases ) )
                return None

            #
            # Count mismatch
            #
            mis_base_U = None
            if int( depth ) > min_depth:

                read_bases = read_bases.replace( '.', ref_base_U )
                read_bases = read_bases.replace( ',', ref_base_U.lower() )

                base_num = {
                    "total_A": 0,
                    "total_C": 0,
                    "total_G": 0,
                    "total_T": 0,
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                }

                #
                # Set data
                #
                data_pair[ data_id ][ 'bases' ] = read_bases
                data_pair[ data_id ][ 'depth' ] = int( depth )

                #
                # Count number
                #
                data_pair[ data_id ][ 'proper_read_depth' ] = 0
                for nuc, qual in zip( read_bases, qual_list ):
                    if not ( qual in filter_quals ):
                        if nuc in 'ATGCatgc':
                            base_num[ nuc ] += 1
                            base_num[ 'total_' + nuc.upper() ] += 1
                            data_pair[ data_id ][ 'proper_read_depth' ] += 1 

                #
                # InsDel
                # Beta distribution
                #
                for type in ( '+', '-' ):
                    if type in indel:
                        for key in indel[ type ].keys():
                            bases = key.split( '\t' )[ 3 ]
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '+' ] = indel[ type ][ key ][ '+' ]
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '-' ] = indel[ type ][ key ][ '-' ]
                            indel_number = \
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 'both' ] = ( indel[ type ][ key ][ '-' ] +
                                                                                           indel[ type ][ key ][ '+' ] )
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '0.1' ] = \
                                scipy.special.btdtri( indel_number + 1, float( depth ) - indel_number + 1, 0.1 )
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 'mid' ] = \
                                ( indel_number + 1 ) / ( float( depth ) + 2 )
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '0.9' ] = \
                                scipy.special.btdtri( indel_number + 1, int( depth ) - indel_number + 1, 0.9 )
                            data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 's_ratio' ] = \
                                float( indel[ type ][ key ][ '+' ] ) / data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 'both' ]

                #
                # skip if reference is 'N'
                #
                if ref_base_U != 'N' and int( data_pair[ data_id ][ 'proper_read_depth' ] ) > min_depth:
                    ref_num = base_num[ 'total_' + ref_base_U ]

                    mis_num = 0
                    for nuc in ( 'A', 'C', 'G', 'T' ):
                        data_pair[ data_id ][ nuc ] = base_num[ nuc ]
                        tmp = nuc.lower()
                        data_pair[ data_id ][ tmp ] = base_num[ tmp ]
                        tmp = 'total_' + nuc
                        data_pair[ data_id ][ tmp ] = base_num[ tmp ]

                        if nuc != ref_base_U:
                            if base_num[ tmp ] > mis_num:
                                mis_num = base_num[ tmp ]
                                mis_base_U = nuc


                    #
                    # Calculate ratio
                    #
                    data_pair[ data_id ][ 'mis_rate' ] = mis_num / float( data_pair[ data_id ][ 'proper_read_depth' ] )
                    data_pair[ data_id ][ 'mis_base' ] = mis_base_U
                    if mis_base_U:
                        data_pair[ data_id ][ 's_ratio' ]  = float( base_num[ mis_base_U ] ) / ( base_num[ mis_base_U ] + base_num[ mis_base_U.lower() ] )
                    else:
                        data_pair[ data_id ][ 's_ratio' ]  = 0

                    #
                    # Beta distribution for SNV
                    #
                    data_pair[ data_id ][ '0.1' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.1 )
                    data_pair[ data_id ][ 'mid' ] = ( mis_num + 1 ) / float( ref_num + mis_num + 2 )
                    data_pair[ data_id ][ '0.9' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.9 )

        #
        # Fisher
        #
        # SNV
        #
        if ( data_pair[ POS_COUNT ] == 2 and
             ref_base_U and
             mis_base_U and
             'mid' in data_pair[ POS_DATA1 ].keys() and
             'mid' in data_pair[ POS_DATA2 ].keys()
           ):
            odds_ratio, fisher_pvalue = fisher(
                        ( ( int( data_pair[ POS_DATA2 ][ 'total_' + ref_base_U ] ),
                            int( data_pair[ POS_DATA1 ][ 'total_' + ref_base_U ] ) ),
                          ( int( data_pair[ POS_DATA2 ][ 'total_' + mis_base_U ] ),
                            int( data_pair[ POS_DATA1 ][ 'total_' + mis_base_U ] ) ) ),
                        alternative='two-sided' )
            data_pair[ POS_FISHER_SNV ] = -math.log( fisher_pvalue, 10 )

        #
        # INDEL
        #
        if data_pair[ POS_COUNT ] == 2 and 'indel' in data_pair[ POS_DATA2 ]:
            fisher_pvalue = None
            for type in data_pair[ POS_DATA2 ][ 'indel' ]:
                for bases in data_pair[ POS_DATA2 ][ 'indel' ][ type ].keys():
                    if type in data_pair[ POS_DATA1 ][ 'indel' ] and bases in data_pair[ POS_DATA1 ][ 'indel' ][ type ]:
                        normal_tmp = data_pair[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] if isinstance( data_pair[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ], int ) else 0 
                        odds_ratio, fisher_pvalue = fisher(
                            ( ( data_pair[ POS_DATA2 ][ 'depth' ] -
                                  data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ],
                                data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ] ),
                              ( data_pair[ POS_DATA1 ][ 'depth' ] - normal_tmp,
                                normal_tmp) ),
                            alternative='two-sided' )

                        if fisher_pvalue != None:
                            if type == '+':
                                data_id = POS_FISHER_INS
                            elif type == '-':
                                data_id = POS_FISHER_DEL

                            if data_pair[ data_id ] == 'N:1.0':
                                data_pair[ data_id ] = bases + ':' + str( -math.log( fisher_pvalue, 10 ) )
                            else:
                                data_pair[ data_id ] += ',' + bases + ':' + str( -math.log( fisher_pvalue, 10 ) )



    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
        data_pair = None
        raise

    return data_pair


############################################################
def print_data( data, w, min_depth, mismatch_rate, posterior_10_quantile, fisher_threshold ):
    str_indel_dict = AutoVivification()
    f_print = False
    f_print_indel = False
    fisher_threshold_log = -math.log( fisher_threshold, 10 )

    outstr = ''

    if data[ POS_COUNT ] == 1:
        #
        # barcode SNV output
        #
        if (
              data[ POS_DATA1 ][ 'mis_rate']      >   mismatch_rate       and
              data[ POS_DATA1 ][ 'mis_base' ]     !=  data[ POS_REF ]     and
              data[ POS_DATA1 ][ 'proper_read_depth' ] >   min_depth      and
              data[ POS_DATA1 ][ '0.1']           >   posterior_10_quantile 
        ):
            f_print = True
            # Genomon output for barcode
            # chr \t start \t end \t ref \t obs \tdepth \t A,C,G,T \t mis \t s_ratio \t 0.1 \t ratio \t 0.9
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6},{7},{8},{9}\t{10:.3f}\t{11:.3f}\t{12:.3f}\t{13:.3f}\t{14:.3f}\n'.format(
                        data[ POS_CHR ],
                        data[ POS_COORD ],
                        data[ POS_COORD ],
                        data[ POS_REF ],
                        data[ POS_DATA1 ][ 'mis_base' ],
                        data[ POS_DATA1 ][ 'depth' ],
                        data[ POS_DATA1 ][ 'total_A' ],
                        data[ POS_DATA1 ][ 'total_C' ],
                        data[ POS_DATA1 ][ 'total_G' ],
                        data[ POS_DATA1 ][ 'total_T' ],
                        data[ POS_DATA1 ][ '0.1' ],
                        data[ POS_DATA1 ][ 'mid' ],
                        data[ POS_DATA1 ][ '0.9' ],
                        data[ POS_DATA1 ][ 'mis_rate' ] if data[ POS_DATA1 ][ 'mis_rate' ] > 0 else '---',
                        data[ POS_DATA1 ][ 's_ratio' ],
                        )

        #
        # InDel output
        #
        if data[ POS_DATA1 ].has_key( 'bases' ):
            indel_number = 0
            for type in data[ POS_DATA1 ][ 'indel' ].keys():
                for bases in data[ POS_DATA1 ][ 'indel' ][ type ].keys():
                    if ( data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] > 0 and
                         data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.1' ]  > posterior_10_quantile ):
                            f_print_indel = True
                            str_indel_dict[ type ][ bases ] = '{0}\t{1}\t{2}\t{3},{4}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\n'.format(
                                        bases if type == '-' else '-',
                                        '-' if type == '-' else bases,
                                        data[ POS_DATA1 ][ 'depth' ],
                                        data[ POS_DATA1 ][ 'depth' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] / float( data[ POS_DATA1 ][ 'depth' ] ),
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '+' ] / float( data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] ),
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.1' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'mid' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.9' ]
                                        )


        if f_print_indel:
            for type in str_indel_dict.keys():
                if type == '+':
                    for bases in str_indel_dict[ type ].keys():
                        #
                        outstr += data[ POS_CHR ] + '\t' +  \
                                  str( data[ POS_COORD ] ) + '\t' + \
                                  str( data[ POS_COORD ] ) + '\t' + \
                                  str_indel_dict[ type ][ bases ]
                elif type == '-':
                    for bases in str_indel_dict[ type ].keys():
                        pos = data[ POS_COORD ]
                        outstr += data[ POS_CHR ] + '\t' + \
                                  str( pos + 1 ) + '\t' + \
                                  str( pos + len( bases ) ) + '\t' + \
                                  str_indel_dict[ type ][ bases ]

    elif data[ POS_COUNT ] == 2:
        #
        # Fisher
        #
        if (
            data[ POS_DATA2 ][ 'mis_rate' ]          > mismatch_rate and
            data[ POS_DATA2 ][ 'proper_read_depth' ] > min_depth     and
            data[ POS_DATA2 ][ 'mis_base' ]    !=  data[ POS_REF ]   and
            data[ POS_FISHER_SNV ]           >  fisher_threshold_log
           ):
            #
            # Genomon output for fisher by comparing nomral and tumor
            # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t A1 \t C1 \t G1 \t T1 \t mis1 \t s_ratio1
            #                        ref2 \t obs2 \tdepth2 \t A2 \t C2 \t G2 \t T2 \t mis2 \t s_ratio2
            #
            f_print = True
            outstr =  '{0}\t{1}\t{2}\t{3}\t{4}'.format(
                        data[ POS_CHR ],
                        data[ POS_COORD ],
                        data[ POS_COORD ],
                        data[ POS_REF ],
                        data[ POS_DATA2 ][ 'mis_base' ],
                        )\
                     + '\t{0},{1},{2},{3}'.format(
                        data[ POS_DATA2 ][ 'total_A' ],
                        data[ POS_DATA2 ][ 'total_C' ],
                        data[ POS_DATA2 ][ 'total_G' ],
                        data[ POS_DATA2 ][ 'total_T' ],
                        )\
                     + '\t{0},{1},{2},{3}'.format(
                        data[ POS_DATA1 ][ 'total_A' ],
                        data[ POS_DATA1 ][ 'total_C' ],
                        data[ POS_DATA1 ][ 'total_G' ],
                        data[ POS_DATA1 ][ 'total_T' ],
                        )\
                     + '\t{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}\n'.format(
                        data[ POS_DATA2 ][ 'mis_rate' ],
                        data[ POS_DATA2 ][ 's_ratio' ],
                        data[ POS_DATA1 ][ 'mis_rate' ],
                        data[ POS_DATA1 ][ 's_ratio' ] if data[ POS_DATA1 ][ 's_ratio' ] > 0 else '---',
                        data[ POS_FISHER_SNV ],
                        )
        else:
            for data_type in ( POS_FISHER_DEL, POS_FISHER_INS ):
                for indel_data in [ x for x in data[ data_type ].split( ',' ) if x != 'N:1.0' ]:
                    bases, fisher_value = indel_data.split( ':' )
                    if float( fisher_value ) > fisher_threshold_log:
                        #
                        # Genomon output for fisher by comparing nomral and tumor
                        # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t obs_depth \t ratio \t 
                        #
                        data_type_symbol = '-' if data_type == POS_FISHER_DEL else '+'
                        disease_misrate = data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ POS_DATA2 ][ 'depth' ] )
                        if disease_misrate > mismatch_rate:
                            f_print = True
                            fisher_tmp_list = indel_data.split( ':' )
                            if isinstance( data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ], int ):
                                normal_tmp = data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] 
                            else:
                                normal_tmp = 0

                            outstr +=  '{0}\t{1}\t{2}\t{3}\t{4}\t'.format(
                                        data[ POS_CHR ],
                                        data[ POS_COORD ] + 1 if data_type == POS_FISHER_DEL else data[ POS_COORD ],
                                        data[ POS_COORD ] + len( bases ) if data_type == POS_FISHER_DEL else data[ POS_COORD ],
                                        fisher_tmp_list[ 0 ] if data_type == POS_FISHER_DEL else '-',
                                        '-' if data_type == POS_FISHER_DEL else fisher_tmp_list[ 0 ]
                                        )\
                                     + '{0},{1}\t{2},{3}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\n'.format(
                                        data[ POS_DATA2 ][ 'depth' ],
                                        data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                        data[ POS_DATA1 ][ 'depth' ],
                                        data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                        disease_misrate,
                                        data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio'],
                                        normal_tmp / float( data[ POS_DATA1 ][ 'depth' ] ),
                                        data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio' ],
                                        fisher_value
                                        )
    if len( outstr ) > 0:
        w.write( outstr )


############################################################
def Pileup_and_count(
        in_bam1 = None,
        in_bam2 = None,
        out_file = None,
        input_mpileup = None,
        ref_fa = None,
        threshold = 15,
        mismatch_rate = 0.07,
        post_10_q = 0.05,
        fisher_threshold = 0.05,
        min_depth = 9,
        compare = False,
        print_header = False
        ):

    global arg
    global target
    global remove_chr
    global filter_quals

    try:

        #
        # Initalize filter quality values
        #
        filter_quals = ''
        for qual in range( 33, 33 + threshold ):
            filter_quals += str( unichr( qual ) )

        #
        # Setup regular expression
        # ([\+\-])[0-9]+[ACGTNacgtn]+
        #
        target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
        remove_chr = re.compile( '\^.' )

        #
        # Open output file and write header
        #
        w = open( out_file, 'w' )

        #
        # Print header only for testing.
        #
        if print_header:
            if compare:
                header_str = "#chr\tstart\tend\tref\tobs\tA,C,G,T\tA,C,G,T\tdis_mis\tdis_s_ration\tctrl_mis\tctrl_s_ratio\tfisher\n"
            else:
                header_str = "#chr\tstart\tend\tref\tobs\tdepth\tA,C,G,T,\tmis\ts_ratio\t0.1\tratio\t0.9\n"

            w.write( header_str )

        #
        # STDIN PIPE
        #   or
        # pysam.mpileup
        #
        if input_mpileup:
            for mpileup in open( arg.input_mpileup, 'rh' ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth, compare )
                if data:
                    print_data( data, w, min_depth, mismatch_rate, post_10_q, fisher_threshold )

        elif in_bam1 and in_bam2:
            for mpileup in pysam.mpileup( '-BQ', '0', '-d', '10000000', '-f', ref_fa, in_bam1, in_bam2 ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth, True )
                if data:
                    print_data( data, w, min_depth, mismatch_rate, post_10_q, fisher_threshold )

        elif in_bam1:
            for mpileup in pysam.mpileup( '-BQ', '0', '-d', '10000000', '-f', ref_fa, in_bam1 ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth, False )
                if data:
                    print_data( data, w, min_depth, mismatch_rate, post_10_q, fisher_threshold )
        else:
            for mpileup in iter( sys.stdin.readline, "" ):
                data = Pileup_out( mpileup, w, threshold, mismatch_rate, min_depth, compare )
                if data:
                    print_data( data, w, min_depth, mismatch_rate, post_10_q, fisher_threshold )

        w.close()

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
        raise

############################################################
def construct_arguments():

    #
    # Arguments
    #
    parser = argparse.ArgumentParser( description = 'Fisher mutation caller',
                                      usage = "%(prog)s [ options ] ")

    parser.add_argument( '-1', '--bam1',          help = '1st bam file ( normal )',  type = str,     default = None )
    parser.add_argument( '-2', '--bam2',          help = '2nd bam file ( disease )', type = str,     default = None )
    parser.add_argument( '-o', '--output',        help = 'Output text file',         type = str,     default = None )

    parser.add_argument( '-i', '--input_mpileup', help = 'Input mpileup file',       type = str,  default = None )
    parser.add_argument( '-c', '--compare',       help = 'Compare two samples', action = 'store_true', default = False )

    parser.add_argument( '-e', '--print_header',  help = 'Print header',        action = 'store_true', default = False )

    parser.add_argument( '-r', '--ref_fa',        help = 'Reference FASTA',          type = str,     default = None )

    parser.add_argument( '-q', '--base_quality',  help = 'Base quality threshold',   type = int,     default = 15 )
    parser.add_argument( '-m', '--mismatch_rate', help = 'Mismatch rate',            type = float,   default = 0.07 )
    parser.add_argument( '-p', '--post_10_q',     help = '10%% posterior quantile threshold', type = float, default = 0.05 )
    parser.add_argument( '-f', '--fisher_value',  help = 'fisher threshold',         type = float,   default = 0.05 )
    parser.add_argument( '-d', '--min_depth',     help = 'Mimimum depth',            type = float,   default = 9 )


    #
    # Log settings
    #
    parser.add_argument( '-g', '--log_file',  help = "Log file name", type = str, default = None )
    parser.add_argument( '-l', '--log_level', help = "Logging level", type = str, default = 'DEBUG' )

    return parser
            

############################################################
def PrintHeader( myself, arg ):
    now = datetime.now()

    logging.info( '#' * 84 )
    logging.info( '# Summary' )
    logging.info( '# Generated by {my}'.format( my = myself ) )
    logging.info( '# %(y)d.%(mo)d.%(d)d.%(h)d:%(mi)d' % { 'y': now.year, 'mo': now.month, 'd': now.day, 'h': now.hour, 'mi': now.minute } )
    logging.info( '#' * 84 + '' )
    logging.info( "bam1: {0}".format( arg.bam1 ) )
    logging.info( "bam2: {0}".format( arg.bam2 ) )
    logging.info( "output: {0}".format( arg.output ) )
    logging.info( "reference fastq: {0}".format( arg.ref_fa ) )
    logging.info( "quality_threshold: {0}".format( arg.base_quality ) )
    logging.info( "mismatch_rate: {0}".format( arg.mismatch_rate ) )
    logging.info( "min_depth: {0}".format( arg.min_depth ) )
    logging.info( '-' * 84 + '' )

############################################################
#
# Main
#
def main():

    #
    # Arguments parse
    #
    argvs = sys.argv
    myself = argvs[ 0 ]
    argc = len(argvs)

    arg_parser = construct_arguments()
    
    if argc < 2:
        arg_parser.print_help()
        sys.exit(1)
                                
    global arg
    arg = arg_parser.parse_args()

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
    level = logging.getLevelName( arg.log_level )

    if arg.log_file:
        logging.basicConfig( filename   = arg.log_file,
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
    Pileup_and_count( 
            in_bam1 = arg.bam1,
            in_bam2 = arg.bam2,
            out_file = arg.output,
            input_mpileup = arg.input_mpileup,
            ref_fa = arg.ref_fa,
            threshold = arg.base_quality,
            mismatch_rate = arg.mismatch_rate,
            post_10_q = arg.post_10_q,
            fisher_threshold = arg.fisher_value,
            min_depth = arg.min_depth,
            compare = arg.compare,
            print_header = arg.print_header
          )

################################################################################
if __name__ == '__main__':
    main()


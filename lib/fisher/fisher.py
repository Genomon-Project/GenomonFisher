#!/usr/bin/python
"""

fisher.py


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
import subprocess

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

def math_log_fisher_pvalue(fisher_pvalue):

    val = float(0.0)
    if fisher_pvalue < 10**(-60):
        val = float(60.0)
    elif fisher_pvalue  > 1.0 - 10**(-10) :
        val = float(0.0)
    else:
        val = -math.log( fisher_pvalue, 10 )
                
    return val


############################################################
def Pileup_out( mpileup, w, min_depth, min_variant_read, compare ):

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
    # if int(mp_list[ 3 ]) < min_depth or ( mp_list_len > 6 and int(mp_list[ 6 ]) < min_depth ):
        return None

    ref_base_plus  = mp_list[ 4 ].count('.')
    ref_base_minus = mp_list[ 4 ].count(',')

    ref_base_count = mp_list[ 4 ].count('.') + mp_list[ 4 ].count(',')
    ins_base_count = mp_list[ 4 ].count('+')
    del_base_count = mp_list[ 4 ].count('-')
    if (int(mp_list[ 3 ]) - ref_base_count + ins_base_count + del_base_count) < min_variant_read:
        return None

    if ref_base_U not in 'ACGTN': return None
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
                  { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'proper_read_depth_plus': 0, 'proper_read_depth_minus': 0, 'proper_read_depth_indel': 0, 'proper_read_depth_indel_plus': 0, 'proper_read_depth_indel_minus': 0,'indel': AutoVivification() },
                  { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'proper_read_depth_plus': 0, 'proper_read_depth_minus': 0, 'proper_read_depth_indel': 0, 'proper_read_depth_indel_plus': 0, 'proper_read_depth_indel_minus': 0,'indel': AutoVivification() },
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
        if int( depth ) >= min_depth:

            read_bases = read_bases.replace( '.', ref_base_U )
            read_bases = read_bases.replace( ',', ref_base_U.lower() )

            base_num = {
                "total_A": 0,
                "total_C": 0,
                "total_G": 0,
                "total_T": 0,
                "total_N": 0,
                "A": 0,
                "C": 0,
                "G": 0,
                "T": 0,
                "N": 0,
                "a": 0,
                "c": 0,
                "g": 0,
                "t": 0,
                "n": 0
            }

            #
            # Set data
            #
            data_pair[ data_id ][ 'bases' ] = read_bases
            data_pair[ data_id ][ 'depth' ] = int( depth )

            #
            # Count number
            #
            for nuc, qual in zip( read_bases, qual_list ):
                if nuc in 'ATGCNacgtn':
                    data_pair[ data_id ][ 'proper_read_depth_indel' ] += 1 
                if nuc in 'ATGCN':
                    data_pair[ data_id ][ 'proper_read_depth_indel_plus' ] += 1 
                if nuc in 'acgtn':
                    data_pair[ data_id ][ 'proper_read_depth_indel_minus' ] += 1 
                if nuc in 'ATGCNacgtn' and not ( qual in filter_quals) :
                    base_num[ nuc ] += 1
                    base_num[ 'total_' + nuc.upper() ] += 1
                if nuc in 'ATGCatgc' and not ( qual in filter_quals):
                    data_pair[ data_id ][ 'proper_read_depth' ] += 1 
                if nuc in 'ATGC' and not ( qual in filter_quals):
                    data_pair[ data_id ][ 'proper_read_depth_plus' ] += 1 
                if nuc in 'atgc' and not ( qual in filter_quals):
                    data_pair[ data_id ][ 'proper_read_depth_minus' ] += 1 

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
                            scipy.special.btdtri( indel_number + 1, float( data_pair[ data_id ][ 'proper_read_depth_indel' ] ) - indel_number + 1, 0.1 )
                        data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 'mid' ] = \
                            ( indel_number + 1 ) / ( float( data_pair[ data_id ][ 'proper_read_depth_indel' ] ) + 2 )
                        data_pair[ data_id ][ 'indel' ][ type ][ bases ][ '0.9' ] = \
                            scipy.special.btdtri( indel_number + 1, int( data_pair[ data_id ][ 'proper_read_depth_indel' ] ) - indel_number + 1, 0.9 )
                        data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 's_ratio' ] = \
                            float( indel[ type ][ key ][ '+' ] ) / data_pair[ data_id ][ 'indel' ][ type ][ bases ][ 'both' ]

            #
            # skip if reference is 'N'
            #
            if ref_base_U != 'N' and int( data_pair[ data_id ][ 'proper_read_depth' ] ) >= min_depth:
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

                if data_id == POS_DATA2 and data_pair[ POS_DATA1 ][ 'mis_base' ]:
                    mis_num = base_num[ 'total_' + data_pair[ POS_DATA1 ][ 'mis_base' ] ]
                    mis_base_U = data_pair[ POS_DATA1 ][ 'mis_base' ]

            ####
                #
                # Calculate ratio
                #
                data_pair[ data_id ][ 'mis_rate' ] = mis_num / float( data_pair[ data_id ][ 'proper_read_depth' ] )
                data_pair[ data_id ][ 'mis_base' ] = mis_base_U
                if mis_base_U and ( base_num[ mis_base_U ] + base_num[ mis_base_U.lower() ] ) > 0:
                    data_pair[ data_id ][ 's_ratio' ]  = float( base_num[ mis_base_U ] ) / ( base_num[ mis_base_U ] + base_num[ mis_base_U.lower() ] )
                # else:
                #    data_pair[ data_id ][ 's_ratio' ]  = float(0)

                #
                # Beta distribution for SNV
                #
                data_pair[ data_id ][ '0.1' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.1 )
                data_pair[ data_id ][ 'mid' ] = ( mis_num + 1 ) / float( ref_num + mis_num + 2 )
                data_pair[ data_id ][ '0.9' ] = scipy.special.btdtri( mis_num + 1, ref_num + 1, 0.9 )

                data_pair[ data_id ][ 'mis_num' ] = mis_num
        
        ###
    #
    # Fisher
    #
    # SNV
    #
    if ( data_pair[ POS_COUNT ] == 2 and
         ref_base_U and
         data_pair[ POS_DATA1 ][ 'mis_base' ] and
         'mid' in data_pair[ POS_DATA1 ].keys() and
         'mid' in data_pair[ POS_DATA2 ].keys() and
         'proper_read_depth' in data_pair[ POS_DATA1 ].keys() and
         'proper_read_depth' in data_pair[ POS_DATA2 ].keys() 
       ):
        odds_ratio, fisher_pvalue = fisher(
                    ( ( int( data_pair[ POS_DATA1 ][ 'total_' + ref_base_U ] ),
                        int( data_pair[ POS_DATA2 ][ 'total_' + ref_base_U ] ) ),
                      ( int( data_pair[ POS_DATA1 ][ 'total_' + data_pair[ POS_DATA1 ][ 'mis_base' ] ] ),
                        int( data_pair[ POS_DATA2 ][ 'total_' + data_pair[ POS_DATA1 ][ 'mis_base' ] ] ) ) ),
                    alternative='two-sided'
                    )

        data_pair[ POS_FISHER_SNV ] = math_log_fisher_pvalue(fisher_pvalue)

    #
    # INDEL
    #
    if ( data_pair[ POS_COUNT ] == 2 and 'indel' in data_pair[ POS_DATA1 ]
       ):
        fisher_pvalue = None
        for type in data_pair[ POS_DATA1 ][ 'indel' ]:
            for bases in data_pair[ POS_DATA1 ][ 'indel' ][ type ].keys():
              
                # if type in data_pair[ POS_DATA1 ][ 'indel' ] and bases in data_pair[ POS_DATA1 ][ 'indel' ][ type ]:

                if not isinstance( data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ], int ):
                    data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ] = 0
                    data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ '+' ] = 0
                    data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ '-' ] = 0

                if (data_pair[ POS_DATA2 ][ 'proper_read_depth_indel' ] >= data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ] and
                    data_pair[ POS_DATA1 ][ 'proper_read_depth_indel' ] >= data_pair[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ]
                    ):
                    odds_ratio, fisher_pvalue = fisher(
                        ( ( data_pair[ POS_DATA1 ][ 'proper_read_depth_indel' ] - data_pair[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ],
                            data_pair[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] ),
                          ( data_pair[ POS_DATA2 ][ 'proper_read_depth_indel' ] - data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ],
                            data_pair[ POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ]) ),
                        alternative='two-sided' )
        
                    if fisher_pvalue != None:
                        if type == '+':
                            data_id = POS_FISHER_INS
                        elif type == '-':
                            data_id = POS_FISHER_DEL

                        if data_pair[ data_id ] == 'N:1.0':
                            data_pair[ data_id ] = bases + ':' + str( math_log_fisher_pvalue(fisher_pvalue) )
                        else:
                            data_pair[ data_id ] += ',' + bases + ':' + str( math_log_fisher_pvalue(fisher_pvalue) )


    return data_pair


############################################################
def print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, posterior_10_quantile, fisher_threshold, min_variant_read ):
    str_indel_dict = AutoVivification()
    f_print_indel = False

    if data[ POS_COUNT ] == 1:
        #
        # barcode SNV output
        #
        if (
              data[ POS_DATA1 ][ 'mis_rate']      >   mismatch_rate_disease and
              data[ POS_DATA1 ][ 'mis_base' ]     !=  data[ POS_REF ]     and
              data[ POS_DATA1 ][ 'proper_read_depth' ] >=  min_depth      and
              data[ POS_DATA1 ][ '0.1']           >   posterior_10_quantile and
              data[ POS_DATA1 ][ 'total_' + data[ POS_DATA1 ][ 'mis_base' ] ] >= min_variant_read
        ):
            # Genomon output for barcode
            # chr \t start \t end \t ref \t obs \tdepth \t A,C,G,T \t mis \t s_ratio \t 0.1 \t ratio \t 0.9
            # outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11:.3f}\t{12:.3f}\t{13:.3f}\t{14:.3f}\t{15:.3f}\n'.format(
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11},{12},{13},{14}\t{15:.3f}\t{16:.3f}\t{17:.3f}\t{18:.3f}\t{19:.3f}\n'.format(
                        data[ POS_CHR ],
                        data[ POS_COORD ],
                        data[ POS_COORD ],
                        data[ POS_REF ],
                        data[ POS_DATA1 ][ 'mis_base' ],
                        data[ POS_DATA1 ][ 'proper_read_depth' ],
                        data[ POS_DATA1 ][ 'total_' + data[ POS_DATA1 ][ 'mis_base' ] ],
                        data[ POS_DATA1 ][ 'proper_read_depth_plus' ],
                        data[ POS_DATA1 ][ data[ POS_DATA1 ][ 'mis_base' ] ],
                        data[ POS_DATA1 ][ 'proper_read_depth_minus' ],
                        data[ POS_DATA1 ][ (data[ POS_DATA1 ][ 'mis_base' ]).lower() ],
                        data[ POS_DATA1 ][ 'total_A' ],
                        data[ POS_DATA1 ][ 'total_C' ],
                        data[ POS_DATA1 ][ 'total_G' ],
                        data[ POS_DATA1 ][ 'total_T' ],
                        data[ POS_DATA1 ][ 'mis_rate' ],
                        data[ POS_DATA1 ][ 's_ratio' ],
                        data[ POS_DATA1 ][ '0.1' ],
                        data[ POS_DATA1 ][ 'mid' ],
                        data[ POS_DATA1 ][ '0.9' ],
                        )
            ###
            w.write( outstr)

        #
        # InDel output
        #
        if data[ POS_DATA1 ].has_key( 'bases' ):
            indel_number = 0
            for type in data[ POS_DATA1 ][ 'indel' ].keys():
                for bases in data[ POS_DATA1 ][ 'indel' ][ type ].keys():
                    if ( data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] > 0 and
                         data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.1' ]  > posterior_10_quantile and
                         data[ POS_DATA1 ][ 'proper_read_depth_indel' ] >=  min_depth and
                        (data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] / float( data[ POS_DATA1 ][ 'depth' ])) > mismatch_rate_disease and
                         data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] >= min_variant_read
                        ):

                            f_print_indel = True
                            # str_indel_dict[ type ][ bases ] = '{0}\t{1}\t{2}\t{3}\t{4},{5}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{10:.3f}\n'.format(
                            str_indel_dict[ type ][ bases ] = '{0}\t{1}\t{2}\t{3}\t{4},{5},{6},{7}\t{8}\t{9:.3f}\t{10:.3f}\t{11:.3f}\t{12:.3f}\t{13:.3f}\n'.format(
                                        bases if type == '-' else '-',
                                        '-' if type == '-' else bases,
                                        data[ POS_DATA1 ][ 'proper_read_depth_indel' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ],
                                        data[ POS_DATA1 ][ 'proper_read_depth_indel_plus' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '+' ],
                                        data[ POS_DATA1 ][ 'proper_read_depth_indel_minus' ],
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ '-' ],
                                        '---',
                                        data[ POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] / float( data[ POS_DATA1 ][ 'proper_read_depth_indel' ] ),
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
                        outstr = data[ POS_CHR ] + '\t' +  \
                                  str( data[ POS_COORD ] ) + '\t' + \
                                  str( data[ POS_COORD ] ) + '\t' + \
                                  str_indel_dict[ type ][ bases ]
                        ###
                        if len( outstr ) > 0:
                            w.write( outstr )

                elif type == '-':
                    for bases in str_indel_dict[ type ].keys():
                        pos = data[ POS_COORD ]
                        #
                        outstr = data[ POS_CHR ] + '\t' + \
                                  str( pos + 1 ) + '\t' + \
                                  str( pos + len( bases ) ) + '\t' + \
                                  str_indel_dict[ type ][ bases ]
                        ###
                        w.write( outstr )

    elif data[ POS_COUNT ] == 2:

        fisher_threshold_log = -math.log( fisher_threshold, 10 )
        #
        # Fisher
        #
        if (
            data[ POS_DATA2 ][ 'mis_rate' ]          < mismatch_rate_normal and
            data[ POS_DATA1 ][ 'mis_rate' ]          > mismatch_rate_disease and
            data[ POS_DATA2 ][ 'proper_read_depth' ] >= min_depth    and
            data[ POS_DATA1 ][ 'proper_read_depth' ] >= min_depth    and
            data[ POS_DATA1 ][ 'mis_base' ]    !=  data[ POS_REF ]   and
            data[ POS_FISHER_SNV ]           >  fisher_threshold_log and
            data[ POS_DATA1 ][ 'total_' + data[ POS_DATA1 ][ 'mis_base' ] ] >= min_variant_read
           ):
            #
            # Genomon output for fisher by comparing nomral and tumor
            # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t A1 \t C1 \t G1 \t T1 \t mis1 \t s_ratio1
            #                        ref2 \t obs2 \tdepth2 \t A2 \t C2 \t G2 \t T2 \t mis2 \t s_ratio2
            #
            ###
            outstr = ('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
                        data[ POS_CHR ],
                        data[ POS_COORD ],
                        data[ POS_COORD ],
                        data[ POS_REF ],
                        data[ POS_DATA1 ][ 'mis_base' ],
                        data[ POS_DATA1 ][ 'proper_read_depth' ],
                        data[ POS_DATA1 ][ 'total_' + data[ POS_DATA1 ][ 'mis_base' ] ],
                        data[ POS_DATA2 ][ 'proper_read_depth' ],
                        data[ POS_DATA2 ][ 'total_' + data[ POS_DATA1 ][ 'mis_base' ] ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ POS_DATA1 ][ 'proper_read_depth_plus' ],
                        data[ POS_DATA1 ][ data[ POS_DATA1 ][ 'mis_base' ] ],
                        data[ POS_DATA1 ][ 'proper_read_depth_minus' ],
                        data[ POS_DATA1 ][ (data[ POS_DATA1 ][ 'mis_base' ]).lower() ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ POS_DATA2 ][ 'proper_read_depth_plus' ],
                        data[ POS_DATA2 ][ data[ POS_DATA1 ][ 'mis_base' ] ],
                        data[ POS_DATA2 ][ 'proper_read_depth_minus' ],
                        data[ POS_DATA2 ][ (data[ POS_DATA1 ][ 'mis_base' ]).lower() ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ POS_DATA1 ][ 'total_A' ],
                        data[ POS_DATA1 ][ 'total_C' ],
                        data[ POS_DATA1 ][ 'total_G' ],
                        data[ POS_DATA1 ][ 'total_T' ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ POS_DATA2 ][ 'total_A' ],
                        data[ POS_DATA2 ][ 'total_C' ],
                        data[ POS_DATA2 ][ 'total_G' ],
                        data[ POS_DATA2 ][ 'total_T' ],
                        )
                     + '\t{0:.3f}'.format(data[ POS_DATA1 ][ 'mis_rate' ])
                     + '\t{0:.3f}'.format(data[ POS_DATA1 ][ 's_ratio' ])
                     + '\t{0:.3f}'.format(data[ POS_DATA2 ][ 'mis_rate' ])
                     )
            if 's_ratio' in data[ POS_DATA2 ].keys():
                outstr += '\t{0:.3f}'.format(data[ POS_DATA2 ][ 's_ratio' ])
            else:
                outstr += '\t---'
            outstr += '\t{0:.3f}'.format(data[ POS_FISHER_SNV ])
                     
            ###
            w.write( outstr +"\n")

        for data_type in ( POS_FISHER_DEL, POS_FISHER_INS ):
            for indel_data in [ x for x in data[ data_type ].split( ',' ) if x != 'N:1.0' ]:
                bases, fisher_value = indel_data.split( ':' )
                if float( fisher_value ) > fisher_threshold_log:
                    #
                    # Genomon output for fisher by comparing nomral and tumor
                    # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t obs_depth \t ratio \t 
                    #
                    data_type_symbol = '-' if data_type == POS_FISHER_DEL else '+'
                    if (data[ POS_DATA2 ][ 'proper_read_depth_indel' ] >= min_depth and
                        data[ POS_DATA1 ][ 'proper_read_depth_indel' ] >= min_depth and
                        (data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ POS_DATA2 ][ 'proper_read_depth_indel' ] )) < mismatch_rate_normal and
                        (data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ POS_DATA1 ][ 'proper_read_depth_indel' ] )) > mismatch_rate_disease and
                        data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] >= min_variant_read
                    ):
                        fisher_tmp_list = indel_data.split( ':' )
                        if isinstance( data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ], int ):
                            normal_tmp = data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] 
                        else:
                            normal_tmp = 0
                        ###
                        outstr = ('{0}\t{1}\t{2}\t{3}\t{4}'.format(
                                    data[ POS_CHR ],
                                    data[ POS_COORD ] + 1 if data_type == POS_FISHER_DEL else data[ POS_COORD ],
                                    data[ POS_COORD ] + len( bases ) if data_type == POS_FISHER_DEL else data[ POS_COORD ],
                                    fisher_tmp_list[ 0 ] if data_type == POS_FISHER_DEL else '-',
                                    '-' if data_type == POS_FISHER_DEL else fisher_tmp_list[ 0 ]
                                    )
                                 + '\t{0}\t{1}\t{2}\t{3}'.format(
                                    data[ POS_DATA1 ][ 'proper_read_depth_indel' ],
                                    data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                    data[ POS_DATA2 ][ 'proper_read_depth_indel' ],
                                    data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                    )
                                 + '\t{0},{1},{2},{3}'.format(
                                    data[ POS_DATA1 ][ 'proper_read_depth_indel_plus' ],
                                    data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ '+' ],
                                    data[ POS_DATA1 ][ 'proper_read_depth_indel_minus' ],
                                    data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ '-' ],
                                    )
                                 + '\t{0},{1},{2},{3}'.format(
                                    data[ POS_DATA2 ][ 'proper_read_depth_indel_plus' ],
                                    data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ '+' ],
                                    data[ POS_DATA2 ][ 'proper_read_depth_indel_minus' ],
                                    data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ '-' ],
                                    )
                                 + '\t{0}\t{1}'.format(
                                    '---',
                                    '---',
                                    )
                                 + '\t{0:.3f}'.format(data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ POS_DATA1 ][ 'proper_read_depth_indel' ] ))
                                 + '\t{0:.3f}'.format(data[ POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio'])
                                 + '\t{0:.3f}'.format(normal_tmp / float( data[ POS_DATA2 ][ 'depth' ] ))
                                 )
                        if 's_ratio' in data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ].keys():
                            outstr += '\t{0:.3f}'.format(data[ POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio' ])
                        else:
                            outstr += '\t---'
                        outstr += '\t{0:.3f}'.format(float(fisher_value))
                        ###
                        w.write( outstr + "\n" )


############################################################
def Pileup_and_count(
        in_bam1,
        in_bam2,
        out_file,
        ref_fa,
        baseq_thres,
        mismatch_rate_disease,
        mismatch_rate_normal,
        post_10_q,
        fisher_threshold,
        min_depth,
        print_header,
        min_variant_read,
        samtools,
        samtools_params,
        region
        ):

    global target
    global remove_chr
    global filter_quals

    #
    # Initalize filter quality values
    #
    filter_quals = ''
    for qual in range( 33, 33 + baseq_thres ):
        filter_quals += str( unichr( qual ) )

    #
    # Setup regular expression
    # ([\+\-])[0-9]+[ACGTNacgtn]+
    #
    target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
    remove_chr = re.compile( '\^.' )

    samtools_params_list = samtools_params.split(" ")

    #
    # Open output file and write header
    #
    w = open( out_file, 'w' )
    FNULL = open(os.devnull, 'w')
    #
    # Print header only for testing.
    #

    if in_bam1 and in_bam2:
        if print_header:
            header_str = "#chr\tstart\tend\tref\talt\tdepth_tumor\tvariantNum_tumor\tdepth_normal\tvariantNum_normal\tbases_tumor\tbases_normal\tA,C,G,T_tumor\tA,C,G,T_normal\tmisRate_tumor\tstrandRatio_tumor\tmisRate_normal\tstrandRatio_normal\tP-value(fisher)\n"
            w.write( header_str )
        cmd_list = [samtools,'mpileup','-f',ref_fa]
        cmd_list.extend(samtools_params_list)
        cmd_list.extend([in_bam1, in_bam2])
        if region:
            cmd_list.insert(2, '-r')
            cmd_list.insert(3, region)
        pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr = FNULL)
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            data = Pileup_out( mpileup, w, min_depth, min_variant_read, True)
            if data:
                print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )

    elif in_bam1 or in_bam2:
        if print_header:
            header_str = "#chr\tstart\tend\tref\talt\tdepth\tvariantNum\tbases\tA,C,G,T\tmisRate\tstrandRatio\t10%_posterior_quantile\tposterior_mean\t90%_posterior_quantile\n"
            w.write( header_str )
        in_bam = in_bam1 if in_bam1 else in_bam2
        cmd_list = [samtools,'mpileup','-f',ref_fa]
        cmd_list.extend(samtools_params_list)
        cmd_list.extend([in_bam])
        if region:
            cmd_list.insert(2, '-r')
            cmd_list.insert(3, region)
        pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr = FNULL)
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            data = Pileup_out( mpileup, w, min_depth, min_variant_read, False )
            if data:
                print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )

    else:
        logging.error( "Input file: {file} not found.".format( file = args.in_bam1 +" "+ args.in_bam2 ) )
        raise
    
    FNULL.close()
    w.close()


#!/usr/bin/python
"""

fisher.py


"""
from __future__ import print_function 
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
from . import print_vcf
from . import print_anno
from . import util
from . import const
import multiprocessing 
import copy
from builtins import chr, str

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
    # mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
    if sys.version_info.major == 3:
        mp_list = mpileup.decode().strip('\n').split( '\t' )
    else:
        mp_list = mpileup.strip('\n').split( '\t' )
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
    # const.POS_CHR = 0
    # const.POS_COORD = 1
    # const.POS_REF = 2
    # const.POS_DATA1 = 3
    # const.POS_DATA2 = 4
    # const.POS_FISHER_SNV = 5
    # const.POS_FISHER_INS = 6
    # const.POS_FISHER_DEL = 7
    # const.POS_COUNT = 8
    #
    data_pair = [ mp_list[ 0 ],
                  int( mp_list[ 1 ] ),
                  mp_list[ 2 ],
                  { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'proper_read_depth_plus': 0, 'proper_read_depth_minus': 0, 'proper_read_depth_indel': 0, 'proper_read_depth_indel_plus': 0, 'proper_read_depth_indel_minus': 0,'indel': util.AutoVivification() },
                  { 'mis_base': ref_base_U, 'mis_rate': 0, 'proper_read_depth': 0, 'proper_read_depth_plus': 0, 'proper_read_depth_minus': 0, 'proper_read_depth_indel': 0, 'proper_read_depth_indel_plus': 0, 'proper_read_depth_indel_minus': 0,'indel': util.AutoVivification() },
                  1.0,
                  'N:1.0',
                  'N:1.0',
                  0 ]


    #
    # Loop for 2 bam file case
    #
    if compare:
        data_pair[ const.POS_COUNT ] = 2
        input_list = [ ( const.POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ),
                       ( const.POS_DATA2, mp_list[ 6 ], mp_list[ 7 ], mp_list[ 8 ] ) ]
    else:
        data_pair[ const.POS_COUNT ] = 1
        input_list = [ ( const.POS_DATA1, mp_list[ 3 ], mp_list[ 4 ], mp_list[ 5 ] ) ]

    #
    # position id,
    # mpileup output 4th row(number of read covering the site),
    # 5th row(read bases),
    # 6th row(base quality)
    #
    for data_id, depth, read_bases, qual_list in input_list:

        indel = util.AutoVivification()

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
        read_bases = read_bases.replace( '$', '' )

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

                if data_id == const.POS_DATA2 and data_pair[ const.POS_DATA1 ][ 'mis_base' ]:
                    mis_num = base_num[ 'total_' + data_pair[ const.POS_DATA1 ][ 'mis_base' ] ]
                    mis_base_U = data_pair[ const.POS_DATA1 ][ 'mis_base' ]

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
    if ( data_pair[ const.POS_COUNT ] == 2 and
         ref_base_U and
         data_pair[ const.POS_DATA1 ][ 'mis_base' ] and
         'mid' in data_pair[ const.POS_DATA1 ].keys() and
         'mid' in data_pair[ const.POS_DATA2 ].keys() and
         'proper_read_depth' in data_pair[ const.POS_DATA1 ].keys() and
         'proper_read_depth' in data_pair[ const.POS_DATA2 ].keys()
       ):
        odds_ratio, fisher_pvalue = fisher(
                    ( ( int( data_pair[ const.POS_DATA1 ][ 'total_' + ref_base_U ] ),
                        int( data_pair[ const.POS_DATA2 ][ 'total_' + ref_base_U ] ) ),
                      ( int( data_pair[ const.POS_DATA1 ][ 'total_' + data_pair[ const.POS_DATA1 ][ 'mis_base' ] ] ),
                        int( data_pair[ const.POS_DATA2 ][ 'total_' + data_pair[ const.POS_DATA1 ][ 'mis_base' ] ] ) ) ),
                    alternative='two-sided'
                    )

        data_pair[ const.POS_FISHER_SNV ] = math_log_fisher_pvalue(fisher_pvalue)

    #
    # INDEL
    #
    if ( data_pair[ const.POS_COUNT ] == 2 and 'indel' in data_pair[ const.POS_DATA1 ]
       ):
        fisher_pvalue = None
        for type in data_pair[ const.POS_DATA1 ][ 'indel' ]:
            for bases in data_pair[ const.POS_DATA1 ][ 'indel' ][ type ].keys():

                # if type in data_pair[ const.POS_DATA1 ][ 'indel' ] and bases in data_pair[ const.POS_DATA1 ][ 'indel' ][ type ]:

                if not isinstance( data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ], int ):
                    data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ] = 0
                    data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ '+' ] = 0
                    data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ '-' ] = 0

                if (data_pair[ const.POS_DATA2 ][ 'proper_read_depth_indel' ] >= data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ] and
                    data_pair[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] >= data_pair[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ]
                    ):
                    odds_ratio, fisher_pvalue = fisher(
                        ( ( data_pair[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] - data_pair[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ],
                            data_pair[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] ),
                          ( data_pair[ const.POS_DATA2 ][ 'proper_read_depth_indel' ] - data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ],
                            data_pair[ const.POS_DATA2 ][ 'indel' ][ type ][ bases ][ 'both' ]) ),
                        alternative='two-sided' )

                    if fisher_pvalue != None:
                        if type == '+':
                            data_id = const.POS_FISHER_INS
                        elif type == '-':
                            data_id = const.POS_FISHER_DEL

                        if data_pair[ data_id ] == 'N:1.0':
                            data_pair[ data_id ] = bases + ':' + str( math_log_fisher_pvalue(fisher_pvalue) )
                        else:
                            data_pair[ data_id ] += ',' + bases + ':' + str( math_log_fisher_pvalue(fisher_pvalue) )


    return data_pair


############################################################
def Pileup_command(
        regions,
        cmd_list,
        min_depth,
        min_variant_read,
        mismatch_rate_disease,
        mismatch_rate_normal,
        post_10_q,
        fisher_threshold,
        is_anno,
        out_file,
        compare_flag,
        w
        ):


    end_idx = 1
    if regions:
        region_list = regions.split(",")   
        end_idx =len(region_list)

    with open(os.devnull, 'w') as FNULL:
    
        for idx in range(end_idx):

            cmd_list_copy = []
            cmd_list_copy = copy.deepcopy(cmd_list)
            if regions:
                cmd_list_copy.insert(2, '-r')
                cmd_list_copy.insert(3, region_list[idx])

            pileup = subprocess.Popen(cmd_list_copy, stdout=subprocess.PIPE)
            for mpileup in pileup.stdout:
                data = Pileup_out( mpileup, w, min_depth, min_variant_read, compare_flag) 
                if data:
                    if is_anno:
                        print_anno.print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )
                    else:
                        print_vcf.print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )
            pileup.stdout.close()  
            pileup.wait()

############################################################
def Pileup_command_multi_thread(
        regions,
        cmd_list,
        min_depth,
        min_variant_read,
        mismatch_rate_disease,
        mismatch_rate_normal,
        post_10_q,
        fisher_threshold,
        is_anno,
        out_file,
        thread_str,
        compare_flag
        ):

    with open(out_file + thread_str, 'w') as w, open(os.devnull, 'w') as FNULL:
        
        region_list = regions.split(",")   
        for idx,target_region in enumerate(region_list):

            cmd_list_copy = []
            cmd_list_copy = copy.deepcopy(cmd_list)
            if target_region:
                cmd_list_copy.insert(2, '-r')
                cmd_list_copy.insert(3, target_region)

            pileup = subprocess.Popen(cmd_list_copy, stdout=subprocess.PIPE)
            for mpileup in pileup.stdout:
                data = Pileup_out( mpileup, w, min_depth, min_variant_read, compare_flag) 
                if data:
                    if is_anno:
                        print_anno.print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )
                    else:
                        print_vcf.print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, min_variant_read )
            pileup.stdout.close()
            pileup.wait()
            
            
############################################################
def Print_header(
        w,
        in_bam1,
        in_bam2,
        sample1,
        sample2,
        ref_fa,
        is_anno
        ):

    #
    # Print metadata only for VCF.
    #
    if not is_anno:
        ref_name, ext = os.path.splitext(ref_fa)
        print_vcf.print_meta(w, ref_name + ".dict", sample1, sample2, in_bam1, in_bam2, ref_fa)

    # 
    # Print header
    #
    if in_bam1 and in_bam2:
        if is_anno:
            print_anno.print_header_pair(w)
        else:
            print_vcf.print_header_pair(w, sample1, sample2)

    elif in_bam1:
        if is_anno:
            print_anno.print_header_single(w)
        else:
            print_vcf.print_header_single(w, sample1)

############################################################
def Pileup_and_count(
        in_bam1,
        in_bam2,
        sample1,
        sample2,
        out_file,
        ref_fa,
        baseq_thres,
        mismatch_rate_disease,
        mismatch_rate_normal,
        post_10_q,
        fisher_threshold,
        min_depth,
        header_flag,
        min_variant_read,
        samtools,
        samtools_params,
        region,
        region_file,
        is_anno
        ):
    global target
    global remove_chr
    global filter_quals

    #
    # Initalize filter quality values
    #
    filter_quals = ''
    for qual in range( 33, 33 + baseq_thres ):
        filter_quals += str( chr( qual ) )

    #
    # Setup regular expression
    # ([\+\-])[0-9]+[ACGTNacgtn]+
    #
    target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
    remove_chr = re.compile( '\^.' )

    samtools_params_list = samtools_params.split(" ")

    region_list = []
    if region_file:
        with open(region_file, 'r') as hin:
            for line in hin:
                region_list.append(line.rstrip('\n'))

    if in_bam1 and in_bam2:
        cmd_list = [samtools,'mpileup','-f',ref_fa]
        cmd_list.extend(samtools_params_list)
        cmd_list.extend([in_bam1, in_bam2])
        compare_flag = True
   
    elif in_bam1 or in_bam2:
        in_bam1 = in_bam1 if in_bam1 else in_bam2
        sample1 = sample1 if sample1 else sample2
        in_bam2 = None
        sample2 = None
        cmd_list = [samtools,'mpileup','-f',ref_fa]
        cmd_list.extend(samtools_params_list)
        cmd_list.extend([in_bam1])
        compare_flag = False

    else:
        logging.error( "Input file: {file} not found.".format( file = in_bam1 +" "+ in_bam2 ) )
        raise ValueError()


    #
    # multi thread
    # 
    if len(region_list) > 0:
        jobs = []
        for idx, target_regions in enumerate(region_list):
            proc = multiprocessing.Process(target = Pileup_command_multi_thread, \
                args = (target_regions, cmd_list, min_depth, min_variant_read, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, is_anno, out_file, "."+str(idx), compare_flag))
            jobs.append(proc)
            proc.start()

        for idx,target_regions in enumerate(region_list):
            jobs[idx].join()

        with open(out_file + ".unsorted", 'w') as w:
            for idx,target_regions in enumerate(region_list):
                with open(out_file +"."+ str(idx), 'r') as hin:
                    for line in hin:
                        print(line.rstrip('\n'), file=w)
         
        subprocess.check_output(["sort", "-k1,1", "-k2,2n", "-V", "-o", out_file+".sorted", out_file+".unsorted"])

        with open(out_file, 'w') as w:
            if header_flag:
                Print_header(w, in_bam1, in_bam2, sample1, sample2, ref_fa, is_anno)
            with open(out_file +".sorted", 'r') as hin:
                for line in hin:
                    print(line.rstrip('\n'), file=w)

        for idx, target_regions in enumerate(region_list):
            os.remove(out_file +"."+ str(idx))
        os.remove(out_file +".sorted")
        os.remove(out_file +".unsorted")
                
    #
    # single thread
    # 
    else:
        with open(out_file, 'w') as w:
            if header_flag:
                Print_header(w, in_bam1, in_bam2, sample1, sample2, ref_fa, is_anno)
            Pileup_command(region, cmd_list, min_depth, min_variant_read, mismatch_rate_disease, mismatch_rate_normal, post_10_q, fisher_threshold, is_anno, out_file, compare_flag, w)
            
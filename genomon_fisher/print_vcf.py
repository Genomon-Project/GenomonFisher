#!/usr/bin/python
"""

print_vcf.py


"""
import sys
import os
import re
import scipy.special
from scipy.stats import fisher_exact as fisher
from . import const
import math

############################################################
const.POS_CHR = 0
const.POS_COORD = 1
const.POS_REF = 2
const.POS_DATA1 = 3
const.POS_DATA2 = 4
const.POS_FISHER_SNV = 5
const.POS_FISHER_INS = 6
const.POS_FISHER_DEL = 7
const.POS_COUNT = 8

############################################################
def print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, posterior_10_quantile, fisher_threshold, min_variant_read ):

    if data[ const.POS_COUNT ] == 1:
        #
        # barcode SNV output
        #
        if (
              data[ const.POS_DATA1 ][ 'mis_rate']      >   mismatch_rate_disease and
              data[ const.POS_DATA1 ][ 'mis_base' ]     !=  data[ const.POS_REF ]     and
              data[ const.POS_DATA1 ][ 'proper_read_depth' ] >=  min_depth      and
              data[ const.POS_DATA1 ][ '0.1']           >   posterior_10_quantile and
              data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ] >= min_variant_read
        ):
            # Genomon output for barcode
            # CHROM{0} \t POS{1} \t ID{2} \t REF{3} \t ALT{4} \t QUAL{5} \t FILTER{6} \t DP={7},AF={9},SB={10} \t FORMAT
            # outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tDP={7};AF={8:.3f};SB={9:.3f}\tDP,DPF,DPR,AD,ADF,ADR,BD1,BDM,BD9\t{7},{10},{11},{12},{13},{14},{15:.3f},{16:.3f},{17:.3f}\n'.format(
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tB10={7:.3f};BM={8:.3f};B90={9:.3f}\tGT:DP:DPF:DPR:AD:ADF:ADR:AF:SB\t{10}:{11}:{12}:{13}:{14}:{15}:{16}:{17:.3f}:{18:.3f}\n'.format(
                        data[ const.POS_CHR ],
                        data[ const.POS_COORD ],
                        '.',
                        data[ const.POS_REF ],
                        data[ const.POS_DATA1 ][ 'mis_base' ],
                        '.',
                        '.',
                        data[ const.POS_DATA1 ][ '0.1' ],
                        data[ const.POS_DATA1 ][ 'mid' ],
                        data[ const.POS_DATA1 ][ '0.9' ],
                        '0/1',
                        data[ const.POS_DATA1 ][ 'proper_read_depth' ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_plus' ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_minus' ],
                        data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ (data[ const.POS_DATA1 ][ 'mis_base' ]).lower() ],
                        data[ const.POS_DATA1 ][ 'mis_rate' ],
                        data[ const.POS_DATA1 ][ 's_ratio' ]
                        )
            ###
            w.write( outstr)

        #
        # InDel output
        #
        if 'bases' in data[ const.POS_DATA1 ]:
            indel_number = 0
            for type in data[ const.POS_DATA1 ][ 'indel' ].keys():
                for bases in data[ const.POS_DATA1 ][ 'indel' ][ type ].keys():
                    if ( data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] > 0 and
                         data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.1' ]  > posterior_10_quantile and
                         data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] >=  min_depth and
                        (data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] / float( data[ const.POS_DATA1 ][ 'depth' ])) > mismatch_rate_disease and
                         data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] >= min_variant_read
                        ):
                            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tB10={7:.3f};BM={8:.3f};B90={9:.3f}\tGT:DP:DPF:DPR:AD:ADF:ADR:AF:SB\t{10}:{11}:{12}:{13}:{14}:{15}:{16}:{17:.3f}:{18:.3f}\n'.format(
                                        data[ const.POS_CHR ],
                                        data[ const.POS_COORD ],
                                        '.',
                                        data[ const.POS_REF ]+bases if type == '-' else data[ const.POS_REF ],
                                        data[ const.POS_REF ] if type == '-' else data[ const.POS_REF ]+bases,
                                        '.',
                                        '.',
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.1' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'mid' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.9' ],
                                        '0/1',
                                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ],
                                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel_plus' ],
                                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel_minus' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '+' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '-' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] / float( data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] ),
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '+' ] / float( data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] )
                                        )
                            ###
                            w.write( outstr)

    elif data[ const.POS_COUNT ] == 2:

        fisher_threshold_log = -math.log( fisher_threshold, 10 )
        #
        # Fisher
        #
        if (
            data[ const.POS_DATA2 ][ 'mis_rate' ]          < mismatch_rate_normal and
            data[ const.POS_DATA1 ][ 'mis_rate' ]          > mismatch_rate_disease and
            data[ const.POS_DATA2 ][ 'proper_read_depth' ] >= min_depth    and
            data[ const.POS_DATA1 ][ 'proper_read_depth' ] >= min_depth    and
            data[ const.POS_DATA1 ][ 'mis_base' ]    !=  data[ const.POS_REF ]   and
            data[ const.POS_FISHER_SNV ]           >  fisher_threshold_log and
            data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ] >= min_variant_read
           ):
            #
            # Genomon output for fisher by comparing nomral and tumor
            # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t A1 \t C1 \t G1 \t T1 \t mis1 \t s_ratio1
            #                        ref2 \t obs2 \tdepth2 \t A2 \t C2 \t G2 \t T2 \t mis2 \t s_ratio2
            # CHROM{0} \t POS{1} \t ID{2} \t REF{3} \t ALT{4} \t QUAL{5} \t FILTER{6} \t DP={7},AF={9},SB={10} \t FORMAT
            #
            ###
            data2_s_ratio = '{0:.3f}'.format(data[ const.POS_DATA2 ][ 's_ratio' ]) if 's_ratio' in data[ const.POS_DATA2 ].keys() else "."
            outstr = ('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tFP={7:.3f}\tGT:DP:DPF:DPR:AD:ADF:ADR:AF:SB'.format(
                        data[ const.POS_CHR ],
                        data[ const.POS_COORD ],
                        '.',
                        data[ const.POS_REF ],
                        data[ const.POS_DATA1 ][ 'mis_base' ],
                        '.',
                        '.',
                        data[ const.POS_FISHER_SNV ]
                        )
                     + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6}:{7:.3f}:{8:.3f}'.format(
                        '0/1',
                        data[ const.POS_DATA1 ][ 'proper_read_depth' ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_plus' ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_minus' ],
                        data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ (data[ const.POS_DATA1 ][ 'mis_base' ]).lower() ],
                        data[ const.POS_DATA1 ][ 'mis_rate' ],
                        data[ const.POS_DATA1 ][ 's_ratio' ]
                        )
                     + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6}:{7:.3f}:{8}'.format(
                        '0/0',
                        data[ const.POS_DATA2 ][ 'proper_read_depth' ],
                        data[ const.POS_DATA2 ][ 'proper_read_depth_plus' ],
                        data[ const.POS_DATA2 ][ 'proper_read_depth_minus' ],
                        data[ const.POS_DATA2 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA2 ][ data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA2 ][ (data[ const.POS_DATA1 ][ 'mis_base' ]).lower() ],
                        data[ const.POS_DATA2 ][ 'mis_rate' ],
                        data2_s_ratio
                        )
                    )

            ###
            w.write( outstr +"\n")

        for data_type in ( const.POS_FISHER_DEL, const.POS_FISHER_INS ):
            for indel_data in [ x for x in data[ data_type ].split( ',' ) if x != 'N:1.0' ]:
                bases, fisher_value = indel_data.split( ':' )
                if float( fisher_value ) > fisher_threshold_log:
                    #
                    # Genomon output for fisher by comparing nomral and tumor
                    # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t obs_depth \t ratio \t
                    #
                    data_type_symbol = '-' if data_type == const.POS_FISHER_DEL else '+'
                    if (data[ const.POS_DATA2 ][ 'proper_read_depth_indel' ] >= min_depth and
                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] >= min_depth and
                        (data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ const.POS_DATA2 ][ 'proper_read_depth_indel' ] )) < mismatch_rate_normal and
                        (data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] )) > mismatch_rate_disease and
                        data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] >= min_variant_read
                    ):
                        fisher_tmp_list = indel_data.split( ':' )
                        if isinstance( data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ], int ):
                            normal_tmp = data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ]
                        else:
                            normal_tmp = 0

                        data2_s_ratio = '{0:.3f}'.format(data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio' ]) if 's_ratio' in data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ].keys() else "."
                        ###
                        outstr = ('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tFP={7:.3f}\tGT:DP:DPF:DPR:AD:ADF:ADR:AF:SB'.format(
                                    data[ const.POS_CHR ],
                                    data[ const.POS_COORD ], 
                                    '.',
                                    data[ const.POS_REF ]+fisher_tmp_list[ 0 ] if data_type == const.POS_FISHER_DEL else data[ const.POS_REF ],
                                    data[ const.POS_REF ] if data_type == const.POS_FISHER_DEL else data[ const.POS_REF ]+fisher_tmp_list[ 0 ],
                                    '.',
                                    '.',
                                    float(fisher_value)
                                    )
                                 + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6}:{7:.3f}:{8:.3f}'.format(
                                    '0/1',
                                    data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ],
                                    data[ const.POS_DATA1 ][ 'proper_read_depth_indel_plus' ],
                                    data[ const.POS_DATA1 ][ 'proper_read_depth_indel_minus' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ '+' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ '-' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ]),
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio']
                                    )
                                 + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6}:{7:.3f}:{8}'.format(
                                    '0/0',
                                    data[ const.POS_DATA2 ][ 'proper_read_depth_indel' ],
                                    data[ const.POS_DATA2 ][ 'proper_read_depth_indel_plus' ],
                                    data[ const.POS_DATA2 ][ 'proper_read_depth_indel_minus' ],
                                    data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                    data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ '+' ],
                                    data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ '-' ],
                                    (normal_tmp / float( data[ const.POS_DATA2 ][ 'depth' ])),
                                    data2_s_ratio,
                                    )
                                )
                        ###
                        w.write( outstr + "\n" )

############################################################
def print_meta(w, ref_dict,sample1,sample2,bam1,bam2,ref_fa):

    # print format
    w.write('##fileformat=VCFv4.2\n')

    # print info and format
    if sample2 != None:
        w.write('##INFO=<ID=FP,Number=1,Type=Float,Description="Minus logarithm of the p-value by Fishers exact test">\n')
    else:
        w.write('##INFO=<ID=B10,Number=1,Type=Float,Description="10% posterior quantile of the beta distribution">\n')
        w.write('##INFO=<ID=BM,Number=1,Type=Float,Description="Posterior mean">\n')
        w.write('##INFO=<ID=B90,Number=1,Type=Float,Description="90% posterior quantile of the beta distribution">\n')
    w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    w.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    w.write('##FORMAT=<ID=DPF,Number=1,Type=Integer,Description="Read depth in the forward strand">\n')
    w.write('##FORMAT=<ID=DPR,Number=1,Type=Integer,Description="Read depth in the reverse strand">\n')
    w.write('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allelic depth">\n')
    w.write('##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Allelic depth in the forward strand">\n')
    w.write('##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Allelic depth in the reverse strand">\n')
    w.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n')
    w.write('##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias">\n')

    # print sample information
    # w.write('##SAMPLE=<ID='+sample1+',SampleName='+sample1+',File='+bam1+'\n')
    # if sample2 != None:
    #     w.write('##SAMPLE=<ID='+sample2+',SampleName='+sample2+',File='+bam2+'\n')

    # print reference information
    with open(ref_dict) as hIN:
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            if F[0] == '@SQ':
                ID = F[1].replace('SN:','')
                length = F[2].replace('LN:','')
                w.write('##contig=<ID='+ID+',length='+length+'>\n')
        w.write('##reference='+ref_fa+'\n')

############################################################
def print_header_pair(w, sample1, sample2):
    w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample1+"\t"+sample2+"\n")

def print_header_single(w, sample):
    w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample+"\n")


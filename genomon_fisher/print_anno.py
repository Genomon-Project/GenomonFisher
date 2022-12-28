#!/usr/bin/python
"""

print_anno.py


"""
import sys
import os
import re
import scipy.special
from scipy.stats import fisher_exact as fisher
import math
from . import util
from . import const

############################################################
def print_data( data, w, min_depth, mismatch_rate_disease, mismatch_rate_normal, posterior_10_quantile, fisher_threshold, min_variant_read, flag_mis_base_0, mismatch_rate_base_0):
    str_indel_dict = util.AutoVivification()
    f_print_indel = False

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
            # chr \t start \t end \t ref \t obs \tdepth \t A,C,G,T \t mis \t s_ratio \t 0.1 \t ratio \t 0.9
            # outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11:.3f}\t{12:.3f}\t{13:.3f}\t{14:.3f}\t{15:.3f}\n'.format(
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7},{8},{9},{10}\t{11},{12},{13},{14}\t{15:.3f}\t{16:.3f}\t{17:.3f}\t{18:.3f}\t{19:.3f}\n'.format(
                        data[ const.POS_CHR ],
                        data[ const.POS_COORD ],
                        data[ const.POS_COORD ],
                        data[ const.POS_REF ],
                        data[ const.POS_DATA1 ][ 'mis_base' ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth' ],
                        data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_plus' ],
                        data[ const.POS_DATA1 ][ data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_minus' ],
                        data[ const.POS_DATA1 ][ (data[ const.POS_DATA1 ][ 'mis_base' ]).lower() ],
                        data[ const.POS_DATA1 ][ 'total_A' ],
                        data[ const.POS_DATA1 ][ 'total_C' ],
                        data[ const.POS_DATA1 ][ 'total_G' ],
                        data[ const.POS_DATA1 ][ 'total_T' ],
                        data[ const.POS_DATA1 ][ 'mis_rate' ],
                        data[ const.POS_DATA1 ][ 's_ratio' ],
                        data[ const.POS_DATA1 ][ '0.1' ],
                        data[ const.POS_DATA1 ][ 'mid' ],
                        data[ const.POS_DATA1 ][ '0.9' ],
                        )
            ###
            w.write( outstr)

        #
        # InDel output
        #
        # print("-----------------------")
        # print(data[ const.POS_DATA1] )
        # print(data)
        # print("-----------------------")
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

                            f_print_indel = True
                            # str_indel_dict[ type ][ bases ] = '{0}\t{1}\t{2}\t{3}\t{4},{5}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{10:.3f}\n'.format(
                            str_indel_dict[ type ][ bases ] = '{0}\t{1}\t{2}\t{3}\t{4},{5},{6},{7}\t{8}\t{9:.3f}\t{10:.3f}\t{11:.3f}\t{12:.3f}\t{13:.3f}\n'.format(
                                        bases if type == '-' else '-',
                                        '-' if type == '-' else bases,
                                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ],
                                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel_plus' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '+' ],
                                        data[ const.POS_DATA1 ][ 'proper_read_depth_indel_minus' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '-' ],
                                        '---',
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] / float( data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] ),
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '+' ] / float( data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'both' ] ),
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.1' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ 'mid' ],
                                        data[ const.POS_DATA1 ][ 'indel' ][ type ][ bases ][ '0.9' ]
                                        )


        if f_print_indel:
            for type in str_indel_dict.keys():
                if type == '+':
                    for bases in str_indel_dict[ type ].keys():
                        #
                        outstr = data[ const.POS_CHR ] + '\t' +  \
                                  str( data[ const.POS_COORD ] ) + '\t' + \
                                  str( data[ const.POS_COORD ] ) + '\t' + \
                                  str_indel_dict[ type ][ bases ]
                        ###
                        if len( outstr ) > 0:
                            w.write( outstr )

                elif type == '-':
                    for bases in str_indel_dict[ type ].keys():
                        pos = data[ const.POS_COORD ]
                        #
                        outstr = data[ const.POS_CHR ] + '\t' + \
                                  str( pos + 1 ) + '\t' + \
                                  str( pos + len( bases ) ) + '\t' + \
                                  str_indel_dict[ type ][ bases ]
                        ###
                        w.write( outstr )

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
            data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ] >= min_variant_read and 
            (data[ const.POS_FISHER_SNV ] >  fisher_threshold_log or 
            (data[ const.POS_DATA2 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ] == 0 and flag_mis_base_0 == True and data[ const.POS_DATA1 ][ 'mis_rate' ] > mismatch_rate_base_0))
           ):
            #
            # Genomon output for fisher by comparing nomral and tumor
            # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t A1 \t C1 \t G1 \t T1 \t mis1 \t s_ratio1
            #                        ref2 \t obs2 \tdepth2 \t A2 \t C2 \t G2 \t T2 \t mis2 \t s_ratio2
            #
            ###
            outstr = ('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
                        data[ const.POS_CHR ],
                        data[ const.POS_COORD ],
                        data[ const.POS_COORD ],
                        data[ const.POS_REF ],
                        data[ const.POS_DATA1 ][ 'mis_base' ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth' ],
                        data[ const.POS_DATA1 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA2 ][ 'proper_read_depth' ],
                        data[ const.POS_DATA2 ][ 'total_' + data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ const.POS_DATA1 ][ 'proper_read_depth_plus' ],
                        data[ const.POS_DATA1 ][ data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA1 ][ 'proper_read_depth_minus' ],
                        data[ const.POS_DATA1 ][ (data[ const.POS_DATA1 ][ 'mis_base' ]).lower() ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ const.POS_DATA2 ][ 'proper_read_depth_plus' ],
                        data[ const.POS_DATA2 ][ data[ const.POS_DATA1 ][ 'mis_base' ] ],
                        data[ const.POS_DATA2 ][ 'proper_read_depth_minus' ],
                        data[ const.POS_DATA2 ][ (data[ const.POS_DATA1 ][ 'mis_base' ]).lower() ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ const.POS_DATA1 ][ 'total_A' ],
                        data[ const.POS_DATA1 ][ 'total_C' ],
                        data[ const.POS_DATA1 ][ 'total_G' ],
                        data[ const.POS_DATA1 ][ 'total_T' ],
                        )
                     + '\t{0},{1},{2},{3}'.format(
                        data[ const.POS_DATA2 ][ 'total_A' ],
                        data[ const.POS_DATA2 ][ 'total_C' ],
                        data[ const.POS_DATA2 ][ 'total_G' ],
                        data[ const.POS_DATA2 ][ 'total_T' ],
                        )
                     + '\t{0:.3f}'.format(data[ const.POS_DATA1 ][ 'mis_rate' ])
                     + '\t{0:.3f}'.format(data[ const.POS_DATA1 ][ 's_ratio' ])
                     + '\t{0:.3f}'.format(data[ const.POS_DATA2 ][ 'mis_rate' ])
                     )
            if 's_ratio' in data[ const.POS_DATA2 ].keys():
                outstr += '\t{0:.3f}'.format(data[ const.POS_DATA2 ][ 's_ratio' ])
            else:
                outstr += '\t---'
            outstr += '\t{0:.3f}'.format(data[ const.POS_FISHER_SNV ])

            ###
            w.write( outstr +"\n")

        for data_type in ( const.POS_FISHER_DEL, const.POS_FISHER_INS ):
            for indel_data in [ x for x in data[ data_type ].split( ',' ) if x != 'N:1.0' ]:
                bases, fisher_value = indel_data.split( ':' )
                data_type_symbol = '-' if data_type == const.POS_FISHER_DEL else '+'
                if (
                    float( fisher_value) >  fisher_threshold_log or 
                    (data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] == 0 and flag_mis_base_0 == True and data[ const.POS_DATA1 ][ 'mis_rate' ] > mismatch_rate_base_0)
                    ):

                    #
                    # Genomon output for fisher by comparing nomral and tumor
                    # chr \t start \t end \t ref1 \t obs1 \tdepth1 \t obs_depth \t ratio \t
                    #
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
                        ###
                        outstr = ('{0}\t{1}\t{2}\t{3}\t{4}'.format(
                                    data[ const.POS_CHR ],
                                    data[ const.POS_COORD ] + 1 if data_type == const.POS_FISHER_DEL else data[ const.POS_COORD ],
                                    data[ const.POS_COORD ] + len( bases ) if data_type == const.POS_FISHER_DEL else data[ const.POS_COORD ],
                                    fisher_tmp_list[ 0 ] if data_type == const.POS_FISHER_DEL else '-',
                                    '-' if data_type == const.POS_FISHER_DEL else fisher_tmp_list[ 0 ]
                                    )
                                 + '\t{0}\t{1}\t{2}\t{3}'.format(
                                    data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                    data[ const.POS_DATA2 ][ 'proper_read_depth_indel' ],
                                    data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ],
                                    )
                                 + '\t{0},{1},{2},{3}'.format(
                                    data[ const.POS_DATA1 ][ 'proper_read_depth_indel_plus' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ '+' ],
                                    data[ const.POS_DATA1 ][ 'proper_read_depth_indel_minus' ],
                                    data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ '-' ],
                                    )
                                 + '\t{0},{1},{2},{3}'.format(
                                    data[ const.POS_DATA2 ][ 'proper_read_depth_indel_plus' ],
                                    data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ '+' ],
                                    data[ const.POS_DATA2 ][ 'proper_read_depth_indel_minus' ],
                                    data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ '-' ],
                                    )
                                 + '\t{0}\t{1}'.format(
                                    '---',
                                    '---',
                                    )
                                 + '\t{0:.3f}'.format(data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 'both' ] / float( data[ const.POS_DATA1 ][ 'proper_read_depth_indel' ] ))
                                 + '\t{0:.3f}'.format(data[ const.POS_DATA1 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio'])
                                 + '\t{0:.3f}'.format(normal_tmp / float( data[ const.POS_DATA2 ][ 'depth' ] ))
                                 )
                        if 's_ratio' in data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ].keys():
                            outstr += '\t{0:.3f}'.format(data[ const.POS_DATA2 ][ 'indel' ][ data_type_symbol ][ bases ][ 's_ratio' ])
                        else:
                            outstr += '\t---'
                        outstr += '\t{0:.3f}'.format(float(fisher_value))
                        ###
                        w.write( outstr + "\n" )


############################################################
def print_header_pair(w):
    w.write("Chr\tStart\tEnd\tRef\tAlt\tdepth_tumor\tvariantNum_tumor\tdepth_normal\tvariantNum_normal\tbases_tumor\tbases_normal\tA_C_G_T_tumor\tA_C_G_T_normal\tmisRate_tumor\tstrandRatio_tumor\tmisRate_normal\tstrandRatio_normal\tP-value(fisher)\n")

def print_header_single(w):
    w.write("Chr\tStart\tEnd\tRef\tAlt\tdepth\tvariantNum\tbases\tA_C_G_T\tmisRate\tstrandRatio\t10%_posterior_quantile\tposterior_mean\t90%_posterior_quantile\n")


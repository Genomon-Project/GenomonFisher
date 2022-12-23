#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .run import run_compare
from .run import run_single

def create_parser():
    prog = "fisher"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
    
    def _create_fisher_pair_parser(subparsers):
        
        fisher_pair_parser = subparsers.add_parser("comparison")
        fisher_pair_parser.add_argument( '-1', '--bam1', help = '1st bam file ( disease )', type = str, default = None, required = True )
        fisher_pair_parser.add_argument( '-2', '--bam2', help = '2nd bam file ( control )', type = str, default = None, required = True )
        fisher_pair_parser.add_argument( '-a', '--sample1', help = '1st sample name ( disease )', type = str, default = None)
        fisher_pair_parser.add_argument( '-b', '--sample2', help = '2nd sample name ( control )', type = str, default = None)
        fisher_pair_parser.add_argument( '-o', '--output', help = 'Output text file', type = str, default = None, required = True)
        fisher_pair_parser.add_argument( '-r', '--ref_fa', help = 'Reference genome in fasta format', type = str, default = None , required = True)
        fisher_pair_parser.add_argument( '-s', '--samtools_path', type = str, default = None, required = True)
        fisher_pair_parser.add_argument( '-S', '--samtools_params', type = str, default = "-q 30 -BQ0 -d 10000000")
        
        fisher_pair_parser.add_argument( '-Q', '--base_quality', help = 'Base quality threshold', type = int, default = 15 )
        fisher_pair_parser.add_argument( '-m', '--min_allele_freq', help = 'minimum amount of disease allele frequency', type = float, default = 0.07 )
        fisher_pair_parser.add_argument( '-M', '--max_allele_freq', help = 'maximum amount of control allele frequency', type = float, default = 0.1 )
        fisher_pair_parser.add_argument( '-f', '--fisher_value', help = 'fisher threshold', type = float, default = 0.05 )
        fisher_pair_parser.add_argument( '-d', '--min_depth', help = 'Mimimum depth', type = float, default = 10 )
        fisher_pair_parser.add_argument( '-v', '--min_variant_read', help = 'Mimimum amount of variant reads (disease)', type = int, default = 4 )
        fisher_pair_parser.add_argument( '-R', '--region', help = 'Region in which pileup is generated', type = str, default = None )
        fisher_pair_parser.add_argument( '-L', '--regions', help = 'The file path of list regions in which pileup is generated', type = str, default = None )
        fisher_pair_parser.add_argument( '-P', '--positions', help = 'The file path of the bed output by pileup.', type = str, default = None )
        fisher_pair_parser.add_argument( '-D', '--flag_mis_base_0', help = 'Ignore the fisher threshold when the number of mismatches is 0 in normal bam.', action = 'store_true', default = False )
        
        fisher_pair_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        fisher_pair_parser.add_argument( '-e', '--print_header', help = 'Print header', action = 'store_true', default = False )
        fisher_pair_parser.add_argument( '-g', '--log_file', help = "Log file name", type = str, default = None )
        fisher_pair_parser.add_argument( '-l', '--log_level', help = "Logging level", type = str, default = 'DEBUG' )
        
        return fisher_pair_parser
        
    def _create_fisher_single_parser(subparsers):
        
        fisher_single_parser = subparsers.add_parser("single")
        fisher_single_parser.add_argument( '-1', '--bam1', help = '1st bam file ', type = str, default = None, required = True )
        fisher_single_parser.add_argument( '-a', '--sample1', help = '1st sample name', type = str, default = None)
        fisher_single_parser.add_argument( '-o', '--output', help = 'Output text file', type = str, default = None, required = True)
        fisher_single_parser.add_argument( '-r', '--ref_fa', help = 'Reference genome in fasta format', type = str, default = None , required = True)
        fisher_single_parser.add_argument( '-s', '--samtools_path', type = str, default = None, required = True)
        fisher_single_parser.add_argument( '-S', '--samtools_params', type = str, default = "-q 30 -BQ0 -d 10000000")
        
        fisher_single_parser.add_argument( '-Q', '--base_quality', help = 'Base quality threshold', type = int, default = 15 )
        fisher_single_parser.add_argument( '-m', '--min_allele_freq', help = 'minimum amount of disease allele frequency', type = float, default = 0.07 )
        fisher_single_parser.add_argument( '-p', '--post_10_q', help = '10 percent posterior quantile threshold', type = float, default = 0.02 )
        fisher_single_parser.add_argument( '-d', '--min_depth', help = 'Mimimum depth', type = float, default = 10 )
        fisher_single_parser.add_argument( '-v', '--min_variant_read', help = 'Mimimum amount of variant reads (disease)', type = int, default = 4 )
        fisher_single_parser.add_argument( '-R', '--region', help = 'region in which pileup is generated', type = str, default = None )
        fisher_single_parser.add_argument( '-L', '--regions', help = 'The file path of list regions in which pileup is generated', type = str, default = None )
        fisher_single_parser.add_argument( '-P', '--positions', help = 'The file path of the bed output by pileup.', type = str, default = None )
        
        fisher_single_parser.add_argument( '-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        fisher_single_parser.add_argument( '-e', '--print_header', help = 'Print header', action = 'store_true', default = False )
        fisher_single_parser.add_argument( '-g', '--log_file', help = "Log file name", type = str, default = None )
        fisher_single_parser.add_argument( '-l', '--log_level', help = "Logging level", type = str, default = 'DEBUG' )


        return fisher_single_parser

    fisher_pair_parser = _create_fisher_pair_parser(subparsers)
    fisher_pair_parser.set_defaults(func = run_compare)
    fisher_single_parser = _create_fisher_single_parser(subparsers)
    fisher_single_parser.set_defaults(func = run_single)
    return parser

#! /usr/bin/env python

import sys
import os
import argparse
import logging
import fisher
import multiprocessing

class ActivePool(object):
    def __init__(self):
        super(ActivePool, self).__init__()
        self.mgr = multiprocessing.Manager()
        self.active = self.mgr.list()
        self.lock = multiprocessing.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)


def fisher_execute(args, target_region, output_file):
        #
        # Main function
        #
        fisher.Pileup_and_count( 
                in_bam1 = args.bam1,
                in_bam2 = args.bam2,
                out_file = output_file,
                ref_fa = args.ref_fa,
                baseq_thres = args.base_quality,
                mismatch_rate_disease = args.min_allele_freq,
                mismatch_rate_normal = args.max_allele_freq,
                post_10_q = None,
                fisher_threshold = args.fisher_value,
                min_depth = args.min_depth,
                min_variant_read = args.min_variant_read,
                samtools = args.samtools_path,
                samtools_params = args.samtools_params,
                region = target_region
              )

def fisher_worker(s, pool, target_region, args):
    name = multiprocessing.current_process().name
    with s:
        pool.makeActive(name)
        fisher_execute(args, target_region, args.output+"."+name)
        pool.makeInactive(name)

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

    header_str = "#chr\tstart\tend\tref\talt\tdepth_tumor\tvariantNum_tumor\tdepth_normal\tvariantNum_normal\tbases_tumor\tbases_normal\tA,C,G,T_tumor\tA,C,G,T_normal\tmisRate_tumor\tstrandRatio_tumor\tmisRate_normal\tstrandRatio_normal\tP-value(fisher)\n"

    pool = ActivePool()
    s = multiprocessing.Semaphore(args.semaphore)
    jobs = []
    if os.path.isfile(args.interval_list):
        i = 1
        with open(args.interval_list, "r") as hin:
            for line in hin:
                jobs.append(multiprocessing.Process(target=fisher_worker, name=str(i), args=(s, pool, line.rstrip("\n"), args)))
                i+=1
    else:
        jobs.append(multiprocessing.Process(target=fisher_worker, name=str(1), args=(s, pool, args.region, args)))

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

    if os.path.isfile(args.interval_list):
        with open(args.output, 'w') as outfile:
            if args.print_header:
                outfile.write(header_str)
            for k in range(1, i):
                with open(args.output+"."+str(k)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.unlink(args.output+"."+str(k))
    else:
        if args.print_header:
            with open(args.output, 'w') as outfile:
                outfile.write(header_str)
                with open(args.output+"."+str(1)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.unlink(args.output+"."+str(1))
        else:
            os.rename(args.output+"."+str(1), args.output)
    

def single_execute(args, target_region, output_file):
    #
    # Main function
    #
    fisher.Pileup_and_count( 
            in_bam1 = args.bam1,
            in_bam2 = None,
            out_file = output_file,
            ref_fa = args.ref_fa,
            baseq_thres = args.base_quality,
            mismatch_rate_disease = args.min_allele_freq,
            mismatch_rate_normal = None,
            post_10_q = args.post_10_q,
            fisher_threshold = None,
            min_depth = args.min_depth,
            min_variant_read = args.min_variant_read,
            samtools = args.samtools_path,
            samtools_params = args.samtools_params,
            region = target_region
          )

def single_worker(s, pool, target_region, args):
    name = multiprocessing.current_process().name
    with s:
        pool.makeActive(name)
        single_execute(args, target_region, args.output+"."+name)
        pool.makeInactive(name)

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

    header_str = "#chr\tstart\tend\tref\talt\tdepth\tvariantNum\tbases\tA,C,G,T\tmisRate\tstrandRatio\t10%_posterior_quantile\tposterior_mean\t90%_posterior_quantile\n"

    pool = ActivePool()
    s = multiprocessing.Semaphore(args.semaphore)
    jobs = []
    if os.path.isfile(args.interval_list):
        i = 1
        with open(args.interval_list, "r") as hin:
            for line in hin:
                jobs.append(multiprocessing.Process(target=single_worker, name=str(i), args=(s, pool, line.rstrip("\n"), args)))
                i+=1
    else:
        jobs.append(multiprocessing.Process(target=single_worker, name=str(1), args=(s, pool, args.region, args)))

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

    if os.path.isfile(args.interval_list):
        with open(args.output, 'w') as outfile:
            if args.print_header:
                outfile.write(header_str)
            for k in range(1, i):
                with open(args.output+"."+str(k)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.unlink(args.output+"."+str(k))
    else:
        if args.print_header:
            with open(args.output, 'w') as outfile:
                outfile.write(header_str)
                with open(args.output+"."+str(1)) as infile:
                    for line in infile:
                        outfile.write(line)
                os.unlink(args.output+"."+str(1))
        else:
            os.rename(args.output+"."+str(1), args.output)



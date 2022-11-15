#! /usr/bin/env python

import sys
import unittest
import os, tempfile, shutil, filecmp
from genomon_fisher import fisher
import genomon_fisher


class TestFisher(unittest.TestCase):

    def setUp(self):
        if sys.version_info.major == 2:
            samtools = "samtools"
        else:
            samtools = "samtools"
        self.samtools = samtools 


    ######################################
    # Tumor/Normal Pair, Annoformat
    ######################################
    def test1(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test1.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.07
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.05
        min_depth = 10
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 30 -BQ0 -d 10000000"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test1.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test2(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test2.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test2.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test3.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test2.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test4.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test4.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test5.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test4.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test6(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test6.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 35
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test6.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test7(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test7.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = "2:10000-243189373"
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test7.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test8(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test8.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = "2:10000-243189373"
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test2.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test9(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test9.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.01
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test9.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test10(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test10.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = cur_dir + "/../data/GRCh37_chr2.bed"
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test7.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    ######################################
    # Tumor Single, Annoformat
    ######################################
    def test11(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test11.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.07
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 10
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 30 -BQ0 -d 10000000"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test11.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test12(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test12.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test12.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test13(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test13.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test13.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test14(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test14.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test14.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test15(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test15.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test15.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test16(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test16.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 35
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test16.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test17(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test17.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = "2:10000-243189373"
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test17.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test18(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test18.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = "2:10000-243189373"
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test13.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test19(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test19.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.1
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test19.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test20(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test20.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = cur_dir + "/../data/GRCh37_chr2.bed"
        is_anno = True
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test17.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    ######################################
    # Tumor/Normal Pair, VCFformat
    ######################################
    def test21(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test21.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.07
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.05
        min_depth = 10
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 30 -BQ0 -d 10000000"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test21.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test22(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test22.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test22.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test23(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test23.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test22.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test24(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test24.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test24.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test25(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test25.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test24.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test26(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test26.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 35
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test26.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test27(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test27.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = "2:10000-243189373"
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test27.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test28(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test28.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.1
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = cur_dir + "/../data/GRCh37_chr2.bed"
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test27.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    ######################################
    # Tumor Single, VCFformat
    ######################################
    def test31(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test31.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.07
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 10
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 30 -BQ0 -d 10000000"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test31.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test32(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test32.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test32.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))
        
    def test33(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test33.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test33.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test34(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test34.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test34.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test35(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test35.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = cur_dir + "/../data/GRCh37_2split.interval_list"
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test35.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test36(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test36.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 35
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = False
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test36.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test37(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test37.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = "2:10000-243189373"
        region_file = None
        positions_bed = None
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test37.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    def test38(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = None
        sample1 = "5929_tumor"
        sample2 = None
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test38.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.02
        mismatch_rate_normal = None
        post_10_q = 0.02
        fisher_threshold = None
        min_depth = 8
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 20 -BQ0 -d 10000000 --ff UNMAP,SECONDARY,QCFAIL,DUP"
        region = None
        region_file = None
        positions_bed = cur_dir + "/../data/GRCh37_chr2.bed"
        is_anno = False
        fisher.Pileup_and_count(in_bam1, in_bam2, sample1, sample2, out_file, ref_fa,
        baseq_thres, mismatch_rate_disease, mismatch_rate_normal, post_10_q,
        fisher_threshold, min_depth, header_flag, min_variant_read, samtools,
        samtools_params, region, region_file, positions_bed, is_anno
        )
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test37.txt"
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

    ######################################
    # Tumor/Normal Pair, Execute
    ######################################
    
    def test41(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        in_bam1 = cur_dir + "/../data/5929_tumor_small.markdup.bam"
        in_bam2 = cur_dir + "/../data/5929_control_small.markdup.bam"
        sample1 = "5929_tumor"
        sample2 = "5929_control"
        out_file =  cur_dir + "/../data/5929_small_mutation_result_test41.txt"
        ref_fa  = cur_dir + "/../data/GRCh37.fa"
        baseq_thres = 15
        mismatch_rate_disease = 0.07
        mismatch_rate_normal = 0.1
        post_10_q = None
        fisher_threshold = 0.05
        min_depth = 10
        header_flag = True
        min_variant_read = 4
        samtools = self.samtools
        samtools_params = "-q 30 -BQ0 -d 10000000"
        region = None
        region_file = None
        positions_bed = None
        is_anno = True
        answer_file = cur_dir + "/../data/5929_small_mutation_result_answer_test1.txt"
        parser = genomon_fisher.parser.create_parser()
        args = parser.parse_args(["comparison", "-1", in_bam1, "-2", in_bam2, "-a", sample1, "-b", sample2, "-o", out_file, "-r", ref_fa, "-s", samtools, "-S", samtools_params, "-Q", str(baseq_thres),  "-m", str(mismatch_rate_disease), "-M", str(mismatch_rate_normal), "-f", str(fisher_threshold), "-d", str(min_depth), "-v", str(min_variant_read), "-e"])
        args.func(args)
        self.assertTrue(filecmp.cmp(out_file, answer_file, shallow=False))

if __name__ == "__main__":
    unittest.main()

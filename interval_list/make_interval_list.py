#!/usr/bin/python

import sys

thread_num = int(sys.argv[1])
gap_file = sys.argv[2]
dict_file = sys.argv[3]
interval_list = sys.argv[4]

genome_dict = {}
with open(dict_file, 'r') as hin:
    for line in hin:
        F = line.split("\t")
        genome_chr = F[1].replace("SN:","")
        genome_len = F[2].replace("LN:","")
        genome_dict[genome_chr] = genome_len
        
tmp_chrom = 0
tmp_start = 0
tmp_end = 0
with open(interval_list, 'w') as hout:
    with open(gap_file, 'r') as hin:
        for line in hin:
            F = line.split("\t")
            chrom = F[0].replace("chr","")
            start = int(F[1])
            end = int(F[2])

            if chrom not in genome_dict: continue

            if chrom == tmp_chrom:
                # if (start+300) > genome_dict[str(chrom)]:
                #    print >> hout, chrom+"\t"+ str(tmp_end-300)+ "\t"+ genome_dict[str(chrom)]
                print >> hout, chrom+"\t"+ str(tmp_end-300)+ "\t"+ str(start+300)

            elif chrom != tmp_chrom:
                if start > 0:
                    print >> hout, chrom+"\t"+ "1" + "\t" + str(start+300)
                if tmp_chrom != 0 and tmp_end != int(genome_dict[str(tmp_chrom)]):
                    print >> hout, tmp_chrom+"\t"+ str(tmp_end-300)+ "\t"+ genome_dict[str(tmp_chrom)]

            tmp_chrom = chrom
            tmp_start = start
            tmp_end = end

pos_list = [''] * thread_num
sum_list = [0] * thread_num
with open(interval_list, 'r') as hin:
    for line in hin:
        F = line.split("\t")
        chrom = F[0]
        start = int(F[1])
        end = int(F[2])
        pos_range = end - start
        min_idx = sum_list.index(min(sum_list))
        pos_list[min_idx] = pos_list[min_idx] +","+chrom+":"+str(start)+"-"+str(end)
        sum_list[min_idx] = sum_list[min_idx] + pos_range

for line in pos_list:
    print line[1:]



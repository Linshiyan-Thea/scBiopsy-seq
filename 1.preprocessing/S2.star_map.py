# ==============================================================================
# Environment: Python 3.12
# No external package dependencies (standard library only: sys, os, re, datetime, argparse)
# External tools: STAR 2.7.11b
# ==============================================================================

import sys
import os
import re
import datetime
import argparse

parser = argparse.ArgumentParser(description='STAR alignment for paired-end RNA-seq reads. \
                                              python S2.star_map_ver0.02.py \
                                              --star_index path/to/STAR/genome/index \
                                              --gtf path/to/annotation.gtf \
                                              Reference: human GRCh38.105 / mouse GRCm38.102')

parser.add_argument('--star_index', type=str, required=True, help='path to the STAR genome index directory (e.g., GRCh38.105 for human, GRCm38.102 for mouse)')
parser.add_argument('--gtf', type=str, required=True, help='path to the GTF annotation file (e.g., Homo_sapiens.GRCh38.105.gtf for human, Mus_musculus.GRCm38.102.gtf for mouse)')
parser.add_argument('-t', '--threads', type=int, default=10, help='number of threads (default: 10)')

args = parser.parse_args()

#define parameters
source_star_path = args.star_index
gtf_file_path = args.gtf
n_threads = args.threads
now_time = datetime.datetime.now()
time_str = now_time.strftime('%Y-%m-%d-%H_%M_%S')
star_out_path = 'star_out_'+time_str
log_out_path = 'log_final_out'
os.system('''mkdir -p ./{path1}/{path2}'''.format(path1 = star_out_path,path2 = log_out_path))
# read files with 1.fq
#os.popeen() will retrun a list of results of shell commands
os.system(' ls *_1_repair_1.fq > clean_file_star_temp.txt')

raw_file = open ('clean_file_star_temp.txt')

for f1 in raw_file:
    f1 = f1.rstrip()
    f2 = re.sub('_1_repair_1.fq','_2_repair_2.fq',f1)
    r1 = f1.split('.')[0]
    r2 = f2.split('.')[0]
    out_prefix = re.sub('_1_repair_1.fq','',f1)
    print('now processing: %s and %s'%(r1,r2))
    os.system('''
STAR --runThreadN {threads} --genomeDir {starpath} --outFileNamePrefix {prefix} --sjdbGTFfile {gtf} --outSAMunmapped Within --readFilesIn {read1}.fq {read2}.fq;
mv {prefix}Log.final.out ./{path1}/{path2};
mv {prefix}Aligned.out.sam {prefix}Log.out {prefix}Log.progress.out {prefix}SJ.out.tab {prefix}_STARgenome ./{path1};
'''.format(threads=n_threads, starpath=source_star_path, read1=r1, read2=r2, prefix=out_prefix, gtf=gtf_file_path, path1=star_out_path, path2=log_out_path)) 
    
os.system('''rm clean_file_star_temp.txt''')

print('star map done!!!')

#!/usr/bin/env python
#coding=utf-8

import os
import re
import argparse

parser = argparse.ArgumentParser(description='Calculate exon and intron mapping rates for all Aligned.out.sam files in the current directory. \
                                              python ExonIntron_RNA_v2.2.py \
                                              -e path/to/exon.bed \
                                              -g path/to/gene.bed')

parser.add_argument('-e', '--exon', type=str, required=True, help='path to the exon BED file (derived from GRCh38.105 for human, GRCm38.102 for mouse)')
parser.add_argument('-g', '--gene', type=str, required=True, help='path to the gene BED file (derived from GRCh38.105 for human, GRCm38.102 for mouse)')

args = parser.parse_args()


def ExonIntron_RNA(sam, exon, gene):
    prefix=re.sub('Aligned.out.sam','',sam)
    os.system("samtools view -@ 15 -bS %s > %s.bam" % (sam, prefix))
    os.system("samtools sort -@ 15 %s.bam -o %s.sorted.bam" % (prefix, prefix))
    os.system("samtools depth %s.sorted.bam > %s.bam.depth" % (prefix, prefix))
    os.system("samtools depth -b %s %s.sorted.bam > %s_exon.bam.depth" % (exon, prefix, prefix))
    os.system("samtools depth -b %s %s.sorted.bam > %s_gene.bam.depth" % (gene, prefix, prefix))

    file_ExonIntron=open(prefix+'_ExonIntron.txt','w')
    file_ExonIntron.write("Sample name: "+prefix+'\n')

    total_bases = 0
    file_total = open(prefix+'.bam.depth','r')
    for line in file_total:
        l = re.sub('\n','',line)
        f = int(l.split('\t')[2])
        total_bases = total_bases + f
    file_ExonIntron.write("Total mapped bases: "+str(total_bases)+'\n')
    file_total.close()
 
    exon_bases = 0
    file_exon = open(prefix+'_exon.bam.depth','r')
    for line in file_exon:
        l = re.sub('\n','',line)
        f = int(l.split('\t')[2])
        exon_bases = exon_bases + f
    p_exon = exon_bases / float(total_bases)
    file_ExonIntron.write("Exon bases: "+str(exon_bases)+'\n')
    file_ExonIntron.write("% of bases mapped to exon: {:.2%}".format(p_exon)+'\n')
    file_exon.close()
    
    gene_bases = 0
    file_gene = open(prefix+'_gene.bam.depth','r')
    for line in file_gene:
        l = re.sub('\n','',line)
        f = int(l.split('\t')[2])
        gene_bases = gene_bases + f
    p_gene = gene_bases / float(total_bases)
    p_intron = p_gene - p_exon
    file_ExonIntron.write("Gene bases: "+str(gene_bases)+'\n')
    file_ExonIntron.write("% of bases mapped to gene: {:.2%}".format(p_gene)+'\n')
    file_ExonIntron.write("% of bases mapped to intron: {:.2%}".format(p_intron)+'\n')
    file_gene.close()
    
    file_ExonIntron.close()
    os.system("rm %s.bam %s.bam.depth %s_exon.bam.depth %s_gene.bam.depth" % (prefix, prefix, prefix, prefix))
    # os.system("rm %s %s.bam %s.bam.depth %s_exon.bam.depth %s_gene.bam.depth" % (sam, prefix, prefix, prefix, prefix))
    print(sam+' is finished!')

exon = args.exon
gene = args.gene
os.system("ls *Aligned.out.sam > list_Aligned_out_sam_inex.tmp.txt")
file_list_sam = open('list_Aligned_out_sam_inex.tmp.txt','r')
for i in file_list_sam:
    sam = re.sub('\n','',i)
    ExonIntron_RNA(sam, exon, gene)
file_list_sam.close()
os.system("rm list_Aligned_out_sam_inex.tmp.txt")

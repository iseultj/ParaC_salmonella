#!/usr/bin/python

# script to merge list of bams with same sample ID - for when you need to process e.g. different library protocols in different ways/ if additional sequencing has been added at a later date. 

from __future__ import division
import sys
import os
import subprocess
from subprocess import call
from subprocess import check_output
#from joblib import Parallel, delayed
import fileinput
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("-p","--parallel",help="Number of samples to run in parallel (default=10)", type=int, default=10)
#REQUIRED NO DEFALUTS
parser.add_argument("-i","--bamlist",help="list of bams",type=str)
parser.add_argument("-o","--outdir",help="output directory",type=str)

#Parse Arguments
args = parser.parse_args()
bamlist  =  open(args.bamlist, 'r')
outdir = str(args.outdir)
print(outdir)
#parallel_jobs = args.parallel

#Get list of samples and bam files
all_files = [line.rstrip() for line in bamlist]
samples = sorted(set([line.split("/")[-1].split('-')[0].split('.')[0].split('_')[0] for line in all_files]))
# function for merging

def merge_duprm(merge_label,file_string):
    merge_command = 'picard MergeSamFiles ' + ' I=' + file_string + ' O=' + outdir + "/" + merge_label + '.mq25.bp34.sorted.grouped.duprm.clipped.bam 2>' + outdir + "/" + merge_label + '.merge.log'
    call(merge_command, shell=True)


for i in samples :
    print(i)
    #Define new sample files from your bam_list and old sample files that might be present in the output directoy
    sample_paths = [file for file in all_files if ("/" + i + "-") in file or ("/" + i + ".") in file or ("/" + i + "_") in file] 
    #Set up list of files for merging with picard
    Merge_Connector = ' I='
    file_string = Merge_Connector.join(sample_paths)
    merge_duprm(i,file_string)

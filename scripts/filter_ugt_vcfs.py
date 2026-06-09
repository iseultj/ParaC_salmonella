#!/usr/bin/python

"""
script for filtering unified genotyper VCFs for min/max DP, maximum heterozygosity and converting to a pseudo-haploid VCF
this is written for pseudogenisation/functional variation relative to a  pangenome reference to try and avoid reference bias
This is designed only for a single-sample VCF, and just prints output to stdout  - so pipe into a bgzip cmd to make a properly compressed VCF.GZ file
Iseult 23.08.2024

Bugfix 17-1-25: make sure DP is int not float

15-12-25 Changed gzip.open(file,'rb') to 'rt' - having problems running on GRACE
"""
# Import statements - iseult remember to add as u go
import argparse
import gzip

parser = argparse.ArgumentParser()

# add arguments
parser.add_argument("-i","--infile", help="input VCF: raw output from unified genotyper with ploidy set to 2")
parser.add_argument("-m","--mean", help="mean coverage from qualimap calculations aligned to linear reference (NOT pangenome ref)")
parser.add_argument("-s","--stdev", help="stdev coverage from qualimap calculations")
parser.add_argument("-c","--minc",help="absolute minimum number of reads supporting a call, default = 3", default = 3)
parser.add_argument("-u","--min_supp",help="minimum proportion of total reads supporting a call, default = 0.9", default = 0.9)


## parse comand line arguments
args = parser.parse_args()

infile = str(args.infile)
meancov = float(args.mean)
stdcov = float(args.stdev)
absmin = float(args.minc)
min_support = float(args.min_supp)
# define minimum , maximum cov 
if absmin > (meancov - 2*stdcov):
    mincov = float(absmin)
else:
    mincov = float(meancov - 2*stdcov)
maxcov = float(meancov + 2*stdcov)

# function to parse variants and filter on coverage and % reads supporting a call; and converting to pseudohaploid.
def variant_filter(variant_line, min_cov, max_cov, min_support):
    """
    Filters each non-header line to remove LowQual, het and sites with coverage outside allowed ranges. Returns pseudohaploid or missing call
    """
    ## parse variant line to get elements of output line
    chrom = variant_line.strip().split()[0]
    pos  = variant_line.strip().split()[1]
    snpid  = variant_line.strip().split()[2]
    ref  = variant_line.strip().split()[3]
    alt  = variant_line.strip().split()[4]
    qual  = variant_line.strip().split()[5]
    filter  = variant_line.strip().split()[6]
    info_string  = variant_line.strip().split()[7] ## might need to come back to this
    fmt_string_list  = variant_line.strip().split()[8].split(":")
    full_gt_list  = variant_line.strip().split()[9].split(":")
    gt_index = fmt_string_list.index("GT")
    
    # 1. if GT = "./."; return line unchanged
    if full_gt_list[gt_index] == "./.":
        output_line = variant_line
    # 2. if filter string is LowQual return empty SNP
    elif filter != ".":
        output_line = chrom + '\t' + pos + '\t' + snpid + '\t' + ref + '\t' + alt + '\t' + "." + '\t' +"." + '\t' +"." + '\t' + "GT" + '\t' + "./."
    else:
        # define DP index
        dp_index = fmt_string_list.index("DP")
        if float(full_gt_list[dp_index]) < float(min_cov) or float(full_gt_list[dp_index]) > float(max_cov):
            output_line = chrom + '\t' + pos + '\t' + snpid + '\t' + ref + '\t' + alt + '\t' + "." + '\t' +"." + '\t' +"." + '\t' + "GT" + '\t' + "./."
        elif full_gt_list[gt_index] == "0/0":
            # variant line remains unchanged
            output_line = variant_line
        elif full_gt_list[gt_index] == "1/1":
            # variant line remains unchanged
            output_line = variant_line
        else:
            # define allele depth filters for het calls
            ad_index =  fmt_string_list.index("AD")
            ad_ref = int(full_gt_list[ad_index].split(',')[0])
            ad_alt  = int(full_gt_list[ad_index].split(',')[1])
            # if ad_ref is > ad_alt and > min AND ad_ref/full_gt_list[dp_index] >= 0.9 set to 0/0
            if ad_ref > ad_alt and float(ad_ref) >= min_cov and  float(ad_ref/float(full_gt_list[dp_index])) >= min_support:
                #info_dict = dict(item.split("=") for item in info_string.split(";"))
                #new_info = 'AN='+ str(info_dict['AN']) + 'DP=' + str(ad_ref) + 'MQ=' + str(info_dict['MQ']) + 'MQ0=' + str(info_dict['MQ0']) #e.g. AN=2;DP=1;MQ=37.00;MQ0=0
                output_line  = chrom + '\t' + pos + '\t' + snpid + '\t' + ref + '\t' + alt + '\t' + qual + '\t' + filter + '\t' + info_string + '\t' + 'GT:DP' + '\t' + '0/0:' + str(ad_ref)
            # elif ad_alt is > ad_ref and > min AND ad_alt/full_gt_list[dp_index] >= 0.9 set to 1/1
            elif ad_alt > ad_ref and float(ad_alt) >= min_cov and float(ad_alt/float(full_gt_list[dp_index])) >= min_support:
                #info_dict = dict(item.split("=") for item in info_string.split(";"))
                #new_info = 'AN='+ str(info_dict['AN']) + 'DP=' + str(ad_alt) + 'MQ=' + str(info_dict['MQ']) + 'MQ0=' + str(info_dict['MQ0'])
                output_line  = chrom + '\t' + pos + '\t' + snpid + '\t' + ref + '\t' + alt + '\t' + qual + '\t' + filter + '\t' + info_string + '\t' + 'GT:DP' + '\t' + '1/1:' + str(ad_alt)
            # else set to missing
            else:
                output_line = chrom + '\t' + pos + '\t' + snpid + '\t' + ref + '\t' + alt + '\t' + "." + '\t' +"." + '\t' +"." + '\t' + "GT" + '\t' + "./."
    return(output_line)



# parse VCF file

with gzip.open(infile, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            # header lines
            print(line.strip())
        else:
            input_line = line.strip()
            output_line = variant_filter(input_line, min_cov=mincov, max_cov=maxcov, min_support=min_support)
            print(output_line)
f.close()

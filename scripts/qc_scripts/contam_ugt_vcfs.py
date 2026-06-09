#!/usr/bin/python

"""
Script for parsing diploid, unfiltered unified genotyper VCFs to assess heteroplasmy in bacterial alignments
"""
# Import statements - iseult remember to add as u go
import argparse
import gzip
import csv
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument("-i","--infile", help="input VCF: raw output from unified genotyper with ploidy set to 2")
parser.add_argument("-o","--outfile",help="prefix for output csv [Default: variant_info]", default="variant_info")


## parse comand line arguments
args = parser.parse_args()

infile = str(args.infile)
outfile = str(args.outfile)

# function to parse variants and filter on coverage and % reads supporting a call; and converting to pseudohaploid.
def variant_parse(variant_line):
    """
    Parses unfiltered VCF to return position, read depth, minor allele support, quality filter in VCF and mutation type
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
    
    # 1. if GT = "./."; skip
    if full_gt_list[gt_index] == "./.":
        outline = [str(pos),"0","0","miss","NA"]
    # 2. if filter string is LowQual return info and flag
    elif filter != ".":
        dp_index = fmt_string_list.index("DP")
        if full_gt_list[gt_index] == "0/0":
            # get DP
            dp = str(full_gt_list[dp_index])
            # minor allele support is 0; no tn/tv info b/c homref
            outline = [str(pos),dp,"0","LowQual","HomRef"]
        elif full_gt_list[gt_index] == "1/1":
            dp = str(full_gt_list[dp_index])
            if (ref == "C" and alt == "T") or (ref=="T" and alt=="C") or (ref == "G" and alt == "A") or (ref == "A" and alt == "G"):
                muttype = "tn"
            else:
                muttype = "tv"
            outline = [str(pos),dp,"0","LowQual",muttype]
        else: # deal with alt calls 
            dp = str(full_gt_list[dp_index])
            ad_index =  fmt_string_list.index("AD")
            ad_ref = int(full_gt_list[ad_index].split(',')[0])
            ad_alt  = int(full_gt_list[ad_index].split(',')[1])
            # if ref > alt, then alt is the minor call
            if ad_ref > ad_alt:
                minor_supp = float(ad_alt)/float(dp)
            # if alt > ref, then ref is the minor call
            else:
                minor_supp = float(ad_ref)/float(dp)
            # define muttype
            if (ref == "C" and alt == "T") or (ref=="T" and alt=="C") or (ref == "G" and alt == "A") or (ref == "A" and alt == "G"):
                muttype = "tn"
            else:
                muttype = "tv"
            outline = [str(pos),dp,str(minor_supp),"LowQual",muttype]
    else: # deal with calls that would actually be used 
        dp_index = fmt_string_list.index("DP")
        if full_gt_list[gt_index] == "0/0":
            # get DP
            dp = str(full_gt_list[dp_index])
            # minor allele support is 0; no tn/tv info b/c homref
            outline = [str(pos),dp,"0",".","HomRef"]
        elif full_gt_list[gt_index] == "1/1":
            dp = str(full_gt_list[dp_index])
            if (ref == "C" and alt == "T") or (ref=="T" and alt=="C") or (ref == "G" and alt == "A") or (ref == "A" and alt == "G"):
                muttype = "tn"
            else:
                muttype = "tv"
            outline = [str(pos),dp,"0",".",muttype]
        else: # deal with alt calls 
            dp = str(full_gt_list[dp_index])
            ad_index =  fmt_string_list.index("AD")
            ad_ref = int(full_gt_list[ad_index].split(',')[0])
            ad_alt  = int(full_gt_list[ad_index].split(',')[1])
            # if ref > alt, then alt is the minor call
            if ad_ref > ad_alt:
                minor_supp = float(ad_alt)/float(dp)
            # if alt > ref, then ref is the minor call
            else:
                minor_supp = float(ad_ref)/float(dp)
            # define muttype
            if (ref == "C" and alt == "T") or (ref=="T" and alt=="C") or (ref == "G" and alt == "A") or (ref == "A" and alt == "G"):
                muttype = "tn"
            else:
                muttype = "tv"
            outline = [str(pos),dp,str(minor_supp),".",muttype]
    return(outline)

 


# parse VCF file
# list of lists for output 
list_of_lists = [["Position","Depth","Minor_Support","QualFilt","MutType"]]
with gzip.open(infile, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            # header lines: skip 
            continue
        else:
            input_line = line.strip()
            output_line = variant_parse(input_line)
            list_of_lists.append(output_line)
f.close()

# write list of lists as CSV file
with open(outfile + ".csv","w") as csvfile:
    stats_writer = csv.writer(csvfile, delimiter = ',')
    stats_writer.writerows(list_of_lists)
    csvfile.close()

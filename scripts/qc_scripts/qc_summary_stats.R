rm(list=ls())
#setwd('~/Desktop/salmonella_scratch/VCF_QC_for_paper')
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

#OptParse parameters
setDTthreads(1)

parser <- OptionParser(description="calculate summary stats and rolling heterozygosity estimates across windows of size X kb (default 10); plots outputs")
parser <- add_option(parser,c("-i","--input"), action="store", help="input CSV file from VCF")
parser <- add_option(parser, c("-s","--size"), action = "store", help="window size for heterozygosity calculations (default 10kb)", default=10000)
parser <- add_option(parser, c("-o","--out"), action = "store", help="output_prefix", default="het_out")
parser <- add_option(parser, c("-n","--step"), action = "store", help="step size for het calculations (default 1kB)", default=1000)
parser <- add_option(parser,c("-m","--min"),action="store",help="absolute minimum coverage to filter on (default = 3)", default = 3)

params <- parse_args(parser)
#### Functions: rolling het calculations ####
rolling_summary <- function(DF, time_col, fun, window_size, step_size, min_window=min(DF$time_col)) {
  # time_col is name of time (position) column
  # fun is function to apply to the subsetted data frames
  # min_window is the start time of the earliest window
  
  times <- DF[, time_col]
  m_t <- max(times) 
  # window_starts is a vector of the windows' minimum times
  window_starts <- seq(from=min_window, to=m_t, by=step_size)
  
  # The i-th element of window_rows is a vector that tells us the row numbers of
  # the data-frame rows that are present in window i 
  window_rows <- lapply(window_starts, function(x) { which(times>=x & times<x+window_size) })
  
  window_summaries <- sapply(window_rows, function(w_r) fun(DF[w_r, ]))
  #print(window_summaries)
  data.frame(start_pos=window_starts, end_pos=window_starts+window_size, summary=window_summaries)
}

het <- function(df) {
  heteroplasmy <- mean(df$Minor_Support, na.rm=T)
  return(heteroplasmy)
}

het_no_md <- function(df) {
  dat_nomd <- df %>% filter(MutType != "tn")
  heteroplasmy <-  mean(dat_nomd$Minor_Support, na.rm=T)
  return(heteroplasmy)
}

#### End Fnx ####

#### Inputs ####
# positions which should be excluded by mask 
to_excl <- fread('~/Desktop/salmonella_scratch/VCF_QC_for_paper/paraC_mask.txt',header=F)
# per individual file
dat <- fread(params$input)
mincov <- as.numeric(params$min) # mimic a normal multiVCF run

#### Initial filters ####

# remove low quality, masked sites and those with quality flagged
dat <- dat %>% filter(!Position %in% to_excl$V1) %>% 
  filter(QualFilt==".") %>% 
  filter(Depth >= mincov)
# work with this data before looking at max coverage filters etc.
#### Summary statistics ####
sumstats <- dat %>% summarise(n_sites=n(),# total sites
                              n_nonref=nrow(dat[MutType!="HomRef"]), # total non-ref
                              n_hetplas=nrow(dat[Minor_Support > 0]), # total het sites
                              n_hetplas_tn=nrow(dat[Minor_Support>0 & MutType=="tn"]),
                              overall_contam=mean(Minor_Support), # this is not weighted by depth - something to be aware off
                              hetpls_contam=mean(Minor_Support[Minor_Support > 0]), # contam estimation from hetpl sites only
                              contam_nomd=mean(Minor_Support[MutType != "tn"]), # contam estimates no md
                              hetpls_contam_nomd=mean(Minor_Support[Minor_Support > 0 & MutType != "tn"]))# contam no md hetpl sites
                  
#### Rolling het calculations ####
df <- data.frame(dat)
rolling_het_with_md <- rolling_summary(df, "Position", fun=het, window_size=params$size, step_size=params$step, min_window=3)
rolling_het_no_md <- rolling_summary(df, "Position", fun=het_no_md, window_size=params$size, step_size=params$step, min_window=3)

colnames(rolling_het_with_md) <- c("start","end","het")
colnames(rolling_het_no_md) <- c("start","end","het_no_md")

merged <- merge(rolling_het_with_md, rolling_het_no_md)
head(merged)
#### Plots ####
# Histogram of minor allele support
hist_plt <- ggplot(dat, aes(x=Minor_Support)) + geom_histogram() + xlim(0.01,0.5) + theme_classic()

rolling_het_plt <- ggplot(merged, aes(x=start,y=het)) + geom_point() + theme_classic()
rolling_het_no_md_plt <- ggplot(merged, aes(x=start,y=het_no_md)) + geom_point() + theme_classic()
#### save plots and dataframes ####
merged$ID <- params$out
sumstats$ID <- params$out
fwrite(sumstats,paste0(params$out,".summary_statistics.csv"))
fwrite(merged,paste0(params$out,".rolling_het.winsize",as.character(params$size),".stepsize",as.character(params$step),".csv"))

pdf(paste0(params$out,"heteroplasmy_histogram_excl_nonhet_sites.pdf"))
hist_plt
dev.off()

pdf(paste0(params$out,".rolling_het.allsites.winsize",as.character(params$size),".stepsize",as.character(params$step),".pdf"))
rolling_het_plt
dev.off()

pdf(paste0(params$out,".rolling_het.nomd_sites.winsize",as.character(params$size),".stepsize",as.character(params$step),".pdf"))
rolling_het_no_md_plt
dev.off()

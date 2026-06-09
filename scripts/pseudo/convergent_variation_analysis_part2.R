rm(list=ls())
wd <- 'convergent_pseudo'
setwd(wd)
# package imports
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggVennDiagram))
#### 0. read in files ####
tr <- read.tree('input_tree_nodeNames.nwk')
# gene-level metadata
meta <- fread('data/pangenome/wgMLST.annotation.csv')
meta <- meta %>% mutate(gene=paste0(`Locus tag`,"_1"))
df <- fread("snppar_annotations.with_fixed_serovar_and_adapted_annotations.csv") # from the part 1 script 
df_internal <- df %>% filter(is_internal==TRUE)
df_internal_max5pc_rev <- df_internal %>% filter(total_rev_nodes/n_descendant_nodes < 0.05)
df_internal_max5pc_rev <- df_internal_max5pc_rev %>% mutate(nonsyn=if_else(impact %in% c("high","missense"),1,0)) %>% 
  mutate(syn=if_else(impact %in% c("synonymous"),1,0))
#### end ####
#### 1. annotation based on ancestral mutations only (not offspring nodes. ) ####
typhi <- c("N4")
# PTA N12
pta_typh <- c("N3")
paratyphiA <- c("N12")
# both paraA and Typhi = N3
# Gallinarum: N40
gallinarum <- c("N40")
#abortusovis N62
abortusovis <- c("N62")
# Fulica N78
fulica <- c("N78")
# abortusequi N86
abortusequi <- c("N86")
choleraesuis_mod <- c("N129")

# choleraesuis with TBS012 on N124
choleraesuis_anc_nodes <- c("N123","N124","N126","N127") 

typhisuis <- c("N140","N139") #shared w knc 139 only 

#modern paratyphi C N157; with Tepos N156 # with TGP N149
# with KSZ N150
paratyphiC_mod <- c("N157")
paratyphiC_with_ancs <- c("N156","N155","N154","N152","N151","N150","N149","N138") # N138 is mrca with typhisuis

infantis <- c("N22")
# Tallahassee N126
tallahassee <- c("N52")
# Bispebjerg N146
bispebjerg <- c("N69")
# Berta N31
berta <- c("N29")
# Birkenhead N42
birkenhead <- c("N97")
## ancients
# need to setdiff with other vectors
# N41: birkenhead with ancients 
bkhd_ancs <- c("N96")
# N52 = bronze age lineage 
ba_ancs <- c("N107")
# N65 KNC clade with XBQM and BPN
knc_cl <- c("N120")
# N51- all ancients except birkenhead, need to setdiff with prev
ancs_except_bkhd <- c("N106","N117","N118","N119","N122")
# ancestral_to_host_spec
ancestral_to_host_spec <- c(ancs_except_bkhd,typhi,paratyphiA,pta_typh,
                            gallinarum,abortusovis,fulica,abortusequi,
                            choleraesuis_mod,choleraesuis_anc_nodes,typhisuis,
                            paratyphiC_mod,paratyphiC_with_ancs)
# ancestral_to_host_generalist
ancestral_to_host_generalist <- c(infantis, tallahassee,bispebjerg,berta,birkenhead)

df_internal_max5pc_rev <- df_internal_max5pc_rev %>% mutate(ancestral_assignment=if_else(node %in% ancestral_to_host_spec, "ancestral_to_host_spec",
                                                                                         if_else(node %in% ancestral_to_host_generalist,"ancestral_to_host_generalist","internal")))

# annotate if synonymous or nonsynonymous
df_internal_max5pc_rev <- df_internal_max5pc_rev %>% 
  mutate(nonsyn=if_else(impact %in% c("high","missense"),1,0)) %>%  
  mutate(syn=if_else(impact %in% c("synonymous"),1,0))
# add var_id and var_id2 for summarising 
df_internal_max5pc_rev <- df_internal_max5pc_rev %>% 
  mutate(varid=paste0(gene,"_", pos, "_", impact))%>% 
  mutate(varid2=paste0(gene,"_", impact))
#### creating summary dataframes ####
# summary of pseudo
gene_level_pseudo_summary <- df_internal_max5pc_rev %>% 
  filter(impact=="high") %>% select(gene, node, serovar, adapted_generalist,ancestral_assignment) %>% 
  unique() %>%
  group_by(gene) %>% summarise(n_pseud_total = n(), 
                               n_pseud_ancestral_hostspec = sum(ancestral_assignment == "ancestral_to_host_spec"),
                               n_pseud_ancestral_generalist = sum(ancestral_assignment=="ancestral_to_host_generalist"),
                               n_pseud_host_specific = sum(adapted_generalist == "host_adapted"), 
                               n_definitely_generalist = sum(serovar %in% c("unassigned_ancestral","infantis","bispebjerg","berta","birkenhead","tallahassee")),
                               clade_labels = paste(serovar, collapse = ";"), 
                               node_labels = paste(node, collapse=";")) %>% 
  mutate(n_ambiguous=n_pseud_total -(n_pseud_host_specific + n_definitely_generalist)) %>% mutate(n_internal=n_pseud_total - (n_ambiguous + n_pseud_ancestral_hostspec + n_pseud_ancestral_generalist))
# potential pseudo recomb list
potential_recomb_pseudo <- df_internal_max5pc_rev %>% filter(impact=="high") %>% select(varid) %>% group_by(varid) %>% summarise(n=n()) %>% filter(n > 1)
# summary excl these sites
gene_level_pseudo_summary_excl_rec <- df_internal_max5pc_rev %>% filter(impact=="high") %>% 
  filter(!varid %in% potential_recomb_pseudo$varid) %>% 
  select(gene, node, serovar, adapted_generalist,ancestral_assignment) %>% 
  unique() %>%
  group_by(gene) %>% summarise(n_pseud_total = n(), 
                               n_pseud_ancestral_hostspec = sum(ancestral_assignment == "ancestral_to_host_spec"),
                               n_pseud_ancestral_generalist = sum(ancestral_assignment=="ancestral_to_host_generalist"),
                               n_pseud_host_specific = sum(adapted_generalist == "host_adapted"), 
                               n_definitely_generalist = sum(serovar %in% c("unassigned_ancestral","infantis","bispebjerg","berta","birkenhead","tallahassee")),
                               clade_labels = paste(serovar, collapse = ";"), 
                               node_labels = paste(node, collapse=";")) %>% 
  mutate(n_ambiguous=n_pseud_total -(n_pseud_host_specific + n_definitely_generalist)) %>% 
  mutate(n_internal=n_pseud_total - (n_ambiguous + n_pseud_ancestral_hostspec + n_pseud_ancestral_generalist))

## same thing for non-synonymous
gene_level_missense_summary <- df_internal_max5pc_rev %>% filter(impact=="missense") %>% select(gene, node, serovar, adapted_generalist,ancestral_assignment) %>% 
  unique() %>%
  group_by(gene) %>% summarise(n_missense_total = n(), 
                               n_missense_ancestral_hostspec = sum(ancestral_assignment == "ancestral_to_host_spec"),
                               n_missense_ancestral_generalist = sum(ancestral_assignment=="ancestral_to_host_generalist"),
                               n_missense_host_specific = sum(adapted_generalist == "host_adapted"), 
                               n_missense_definitely_generalist = sum(serovar %in% c("unassigned_ancestral","infantis","bispebjerg","berta","birkenhead","tallahassee")),
                               clade_labels_missense = paste(serovar, collapse = ";"), node_labels_missense = paste(node, collapse=";")) %>% 
  mutate(n_ambiguous_missense=n_missense_total -(n_missense_host_specific + n_missense_definitely_generalist)) %>% 
  mutate(n_internal_missense=n_missense_total - (n_ambiguous_missense + n_missense_ancestral_hostspec + n_missense_ancestral_generalist))
# potential pseudo recomb list
potential_recomb_missense <- df_internal_max5pc_rev %>% filter(impact=="missense") %>% select(varid) %>% group_by(varid) %>% summarise(n=n()) %>% filter(n > 1)
# summary excl these sites
gene_level_missense_summary_excl_rec <- df_internal_max5pc_rev %>% filter(impact=="missense") %>% 
  filter(!varid %in% potential_recomb_missense$varid) %>% 
  select(gene, node, serovar, adapted_generalist,ancestral_assignment) %>% 
  unique() %>%
  group_by(gene) %>% summarise(n_missense_total = n(), 
                               n_missense_ancestral_hostspec = sum(ancestral_assignment == "ancestral_to_host_spec"),
                               n_missense_ancestral_generalist = sum(ancestral_assignment=="ancestral_to_host_generalist"),
                               n_missense_host_specific = sum(adapted_generalist == "host_adapted"), 
                               n_missense_definitely_generalist = sum(serovar %in% c("unassigned_ancestral","infantis","bispebjerg","berta","birkenhead","tallahassee")),
                               clade_labels_missense = paste(serovar, collapse = ";"), node_labels_missense = paste(node, collapse=";")) %>% 
  mutate(n_ambiguous_missense=n_missense_total -(n_missense_host_specific + n_missense_definitely_generalist)) %>% 
  mutate(n_internal_missense=n_missense_total - (n_ambiguous_missense + n_missense_ancestral_hostspec + n_missense_ancestral_generalist))
## same thing for synonymous

gene_level_syn_summary <- df_internal_max5pc_rev %>% filter(impact=="synonymous") %>% select(gene, node, serovar, adapted_generalist,ancestral_assignment) %>% 
  unique() %>%
  group_by(gene) %>% summarise(n_syn_total = n(), n_syn_host_specific = sum(adapted_generalist == "host_adapted"), 
                               n_syn_ancestral_hostspec = sum(ancestral_assignment == "ancestral_to_host_spec"),
                               n_syn_ancestral_generalist = sum(ancestral_assignment=="ancestral_to_host_generalist"),
                               n_syn_definitely_generalist = sum(serovar %in% c("unassigned_ancestral","infantis","bispebjerg","berta","birkenhead","tallahassee")),
                               clade_labels_syn = paste(serovar, collapse = ";"), node_labels_syn = paste(node, collapse=";")) %>% 
  mutate(n_ambiguous_syn=n_syn_total -(n_syn_host_specific + n_syn_definitely_generalist)) %>%
  mutate(n_internal_syn=n_syn_total - (n_ambiguous_syn + n_syn_ancestral_hostspec + n_syn_ancestral_generalist))
# potential pseudo recomb list
potential_recomb_syn <- df_internal_max5pc_rev %>% filter(impact=="synonymous") %>% select(varid) %>% group_by(varid) %>% summarise(n=n()) %>% filter(n > 1)
# summary excl these sites
gene_level_syn_summary_excl_rec <- df_internal_max5pc_rev %>% filter(impact=="synonymous") %>% 
  filter(!varid %in% potential_recomb_syn$varid) %>% 
  select(gene, node, serovar, adapted_generalist,ancestral_assignment) %>% 
  unique() %>%
  group_by(gene) %>% summarise(n_syn_total = n(), n_syn_host_specific = sum(adapted_generalist == "host_adapted"), 
                               n_syn_ancestral_hostspec = sum(ancestral_assignment == "ancestral_to_host_spec"),
                               n_syn_ancestral_generalist = sum(ancestral_assignment=="ancestral_to_host_generalist"),
                               n_syn_definitely_generalist = sum(serovar %in% c("unassigned_ancestral","infantis","bispebjerg","berta","birkenhead","tallahassee")),
                               clade_labels_syn = paste(serovar, collapse = ";"), node_labels_syn = paste(node, collapse=";")) %>% 
  mutate(n_ambiguous_syn=n_syn_total -(n_syn_host_specific + n_syn_definitely_generalist)) %>%
  mutate(n_internal_syn=n_syn_total - (n_ambiguous_syn + n_syn_ancestral_hostspec + n_syn_ancestral_generalist))
#### end summary dfs ####

# summary dataframes to join
# with potential recombination
join_all <- purrr::reduce(list(gene_level_pseudo_summary, gene_level_missense_summary, gene_level_syn_summary), full_join, by="gene")

# without potential recombination
join_excl_recomb <- purrr::reduce(list(gene_level_pseudo_summary_excl_rec,gene_level_missense_summary_excl_rec, gene_level_syn_summary_excl_rec),full_join,by="gene")

# merge with gene level metadata
meta <- fread('data/pangenome/wgMLST.annotation.csv')
meta <- meta %>% mutate(gene=paste0(`Locus tag`,"_1"))
join_all <- merge(join_all,meta,by="gene",all.x=T)
join_excl_recomb <- merge(join_excl_recomb,meta,by="gene",all.x=T)
# set NA to zero
join_all[is.na(join_all)] <- 0
join_excl_recomb[is.na(join_excl_recomb)] <- 0
# save these files 
# fwrite(join_all,"gene_level_summaries.all_mut_types.w_possible_recomb.annotated_correctly.with_gene_meta.tsv")
# fwrite(join_excl_recomb,"gene_level_summaries.all_mut_types.recombfilt.annotated_correctly.with_gene_meta.tsv")
join_all <- fread('gene_level_summaries.all_mut_types.w_possible_recomb.annotated_correctly.with_gene_meta.tsv')
join_excl_recomb <- fread("gene_level_summaries.all_mut_types.recombfilt.annotated_correctly.with_gene_meta.tsv")
#### normalisation factors ####
#sum_pseudo/sum_syn or sum_missense/sum_syn

norm.pseudo.rec <- sum(join_all$n_pseud_total, na.rm=T)/sum(join_all$n_syn_total, na.rm=T)
norm.pseudo.recfilt <- sum(join_excl_recomb$n_pseud_total, na.rm=T)/sum(join_excl_recomb$n_syn_total, na.rm=T)
norm.missense.rec <- sum(join_all$n_missense_total, na.rm=T)/sum(join_all$n_syn_total, na.rm=T)
norm.missense.recfilt <- sum(join_excl_recomb$n_missense_total, na.rm=T)/sum(join_excl_recomb$n_syn_total, na.rm=T)


## plotting obs vs exp
ggplot(join_all, aes(x=n_syn_total*norm.pseudo.rec ,y=n_pseud_total,fill=n_pseud_host_specific,col=n_pseud_host_specific)) + 
  geom_point() + theme_bw() 
## look at this in context of ancestral only
ggplot(join_all, aes(x=n_syn_total*norm.pseudo.rec ,y=n_pseud_total,
                     fill=n_pseud_ancestral_hostspec,col=n_pseud_ancestral_hostspec)) + 
  geom_point() + theme_bw()

# excl recomb
ggplot(join_excl_recomb, aes(x=n_syn_total*norm.pseudo.recfilt ,
                             y=n_pseud_total,
                             fill=n_pseud_host_specific,col=n_pseud_host_specific)) + 
  geom_point() + theme_bw() 

ggplot(join_excl_recomb, aes(x=n_syn_total*norm.pseudo.recfilt ,
                             y=n_pseud_total,
                             fill=n_pseud_ancestral_hostspec,col=n_pseud_ancestral_hostspec)) + 
  geom_point() + theme_bw() 
# missense

ggplot(join_all, aes(x=n_syn_total*norm.missense.rec ,
                     y=n_missense_total,
                     fill=n_missense_host_specific,col=n_missense_host_specific)) + 
  geom_point() + theme_bw() 
ggplot(join_all, aes(x=n_syn_total*norm.missense.rec ,
                     y=n_missense_total,
                     fill=n_missense_ancestral_hostspec,col=n_missense_ancestral_hostspec)) + 
  geom_point() + theme_bw() 

# excl recomb
ggplot(join_excl_recomb, aes(x=n_syn_total*norm.missense.recfilt 
                             ,y=n_missense_total,
                             fill=n_missense_host_specific,
                             col=n_missense_host_specific)) + 
  geom_point() + theme_bw()

ggplot(join_excl_recomb, aes(x=n_syn_total*norm.missense.recfilt 
                             ,y=n_missense_total,
                             fill=n_missense_ancestral_hostspec,
                             col=n_missense_ancestral_hostspec)) + 
  geom_point() + theme_bw() 
# missense looks like a reasonable 1-1 relationship; harder to see with pseudo

# getting gene-level distribution of deviation from expectations
# step 1 - set NA to 0 for counts where necessary
join_all <- join_all %>% map_if(is.numeric,~ifelse(is.na(.x),0,.x)) %>% as.data.table()
# similarly for those without potential recombination
join_excl_recomb <- join_excl_recomb  %>% map_if(is.numeric,~ifelse(is.na(.x),0,.x)) %>% as.data.table()
ggplot(join_all, aes(x=n_pseud_total - (n_syn_total*norm.pseudo.rec))) + geom_density() + theme_bw()
# a bit bumpy; long tail
ggplot(join_all, aes(x=n_missense_total - (n_syn_total*norm.missense.rec))) + geom_density() + theme_bw()


# looks more like a normal distribution for missense, but again a bit of a tail indicating excess missense in a small number of genes

ggplot(join_excl_recomb, aes(x=n_pseud_total - (n_syn_total*norm.pseudo.recfilt))) + 
  geom_density() + theme_bw()
# a bit bumpier because fewer data points
ggplot(join_excl_recomb, aes(x=n_missense_total - (n_syn_total*norm.missense.recfilt))) + 
  geom_density() + theme_bw()

# get quantiles for filtering based on total pseudo/syn counts
quantile((join_all$n_pseud_total - join_all$n_syn_total*norm.pseudo.rec), c(0.025,0.975)) # -0.460714  1.756093 
quantile((join_all$n_missense_total - join_all$n_syn_total*norm.missense.rec), c(0.025,0.975)) #-5.834255  7.570356 

quantile((join_excl_recomb$n_pseud_total - join_excl_recomb$n_syn_total*norm.pseudo.recfilt), c(0.025,0.975)) #-0.4783817  1.6412137
quantile((join_excl_recomb$n_missense_total - join_excl_recomb$n_syn_total*norm.missense.recfilt), c(0.025,0.975)) # -5.253978  5.552020
# then filter for 0 definitely host generalist - refined list of candidates.
## add a column to your df for delta_pseudo and delta_missense

join_all <- join_all %>% mutate(delta_pseudo = (n_pseud_total - (n_syn_total*norm.pseudo.rec)),
                                delta_missense = (n_missense_total - (n_syn_total*norm.missense.rec)))
join_excl_recomb <- join_excl_recomb %>% mutate(delta_pseudo = (n_pseud_total - (n_syn_total*norm.pseudo.recfilt)),
                                                delta_missense = (n_missense_total - (n_syn_total*norm.missense.recfilt)))

# filtering for top 2.5% pseudo only (i don't think the nonsynonymous metric makes a huge amount of sense when not considering identical mutations)
# it's more of a sanity check
# then, with top 2.5% pseudo, get all with 0 known host generalists
top_pseudo_hostspec_norecombfilt <- join_all %>% filter(delta_pseudo >1.756093  ) %>% 
  filter(n_definitely_generalist==0)
top_pseudo_hostspec_recombfilt <- join_excl_recomb %>% filter(delta_pseudo >1.6412137) %>% 
  filter(n_definitely_generalist==0)
nrow(top_pseudo_hostspec_norecombfilt) # 47
nrow(top_pseudo_hostspec_recombfilt) # 58

# save
fwrite(top_pseudo_hostspec_norecombfilt,"top_pseudo_hostspec_norecombfilt.tsv",sep="\t")
fwrite(top_pseudo_hostspec_recombfilt,"top_pseudo_hostspec_recombfilt.tsv",sep="\t")
## also save all host specific pseudo regardless of delta_pseudo
all_hostspec_pseudo_with_rec <- join_all %>% filter(n_pseud_total > 1) %>% filter(n_definitely_generalist==0)
all_hostspec_pseudo_recfilt <- join_excl_recomb %>% filter(n_pseud_total > 1)%>% filter(n_definitely_generalist==0)

nrow(all_hostspec_pseudo_recfilt) # 87
nrow(all_hostspec_pseudo_with_rec) # 104
# fwrite(all_hostspec_pseudo_with_rec,"all_pseudo_hostspec_norecombfilt.tsv",sep="\t")
# fwrite(all_hostspec_pseudo_recfilt,"all_pseudo_hostspec_recombfilt.tsv",sep="\t")
#### end of this section #### 

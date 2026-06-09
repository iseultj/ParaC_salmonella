## script for post-processing of SNPPAR results, once variants are annotated based on revertant status and % offspring which are revertant
rm(list=ls())
wd <- 'convergent_pseudo'
setwd(wd)
# package imports
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggVennDiagram))
#### 0. setup ####
# annotated tree 
tr <- read.tree('input_tree_nodeNames.nwk') # from snppar output
# gene-level metadata
meta <- fread('data/pangenome/wgMLST.annotation.csv')
meta <- meta %>% mutate(gene=paste0(`Locus tag`,"_1"))
# SNPPAR results
high <- fread('full_info_from_snppar.high_impact.annos.txt') # from snppar output
head(high) #gene   pos   node is_internal is_fixed n_descendant_nodes revcount revertant_nodes total_rev_nodes
high$impact <- "high"
miss <- fread("full_info_from_snppar.missense.annos.txt") # from snppar output
miss$impact <- "missense"
syn <- fread("full_info_from_snppar.syn.annos.txt") # from snppar output
syn$impact <- "synonymous"

df <- rbindlist(list(high,miss,syn))
#### end setup ####  
#### 1. gene-level nonsynonymous to synonymous ratios ####
# get degeneracy data
degeneracy_data <- fread('data/pangenome/possible_syn_nonsyn_events.by_gene.csv')
# filter df to restrict to internal variation
internal <- df %>% filter(is_internal==TRUE)
unique(internal)
internal <- internal %>% mutate(nonsyn=if_else(impact %in% c("high","missense"),1,0)) %>% 
  mutate(syn=if_else(impact %in% c("synonymous"),1,0))
#turn internal into gene + syn count/ non syn count
internal <- internal %>% group_by(gene) %>% mutate(nonsyn_count=sum(nonsyn), syn_count=sum(syn)) %>% 
  select(gene, nonsyn_count, syn_count) %>% unique()
# add degeneracy data information
colnames(degeneracy_data) <- c("gene","possible_syn","possible_nonsyn")
internal <- merge(internal, degeneracy_data)
internal <- internal %>% mutate(dN=(nonsyn_count/possible_nonsyn), dS=(syn_count/possible_syn)) %>% mutate(dNdS=(dN/dS))
ggplot(internal, aes(x=dN, y=dS, fill=dNdS, col=dNdS)) + 
  geom_point() + theme_bw()
labels <- internal %>% filter(dNdS >1)
labels <- labels[!is.infinite(labels$dNdS),]

ggplot(internal, aes(y=dN, x=dS, fill=dNdS, col=dNdS)) + 
  geom_point() + theme_bw() + 
  geom_label(data=labels, aes(label=gene), 
             col="black",fill=alpha("white",0.2))
internal_2 <- internal %>% filter(dS > 0)
quantile(internal_2$dNdS, na.rm=T, c(0.9,0.99,0.999))
pc99 <- internal_2 %>% filter(dNdS > 1.2941176)
tmp1 <- meta %>% filter(gene %in% pc99$gene)
tmp1 # mainly hypothetical proteins, CheW , membrane proteins etc. 
# internal and fixed
internal_fixed <- df %>% filter(is_internal==TRUE)%>% filter(is_fixed==TRUE)
internal_fixed <- internal_fixed %>% mutate(nonsyn=if_else(impact %in% c("high","missense"),1,0)) %>% 
  mutate(syn=if_else(impact %in% c("synonymous"),1,0))
#turn internal into gene + syn count/ non syn count
internal_fixed <- internal_fixed %>% group_by(gene) %>% mutate(nonsyn_count=sum(nonsyn), syn_count=sum(syn)) %>% 
  select(gene, nonsyn_count, syn_count) %>% unique()
# add degeneracy data information

internal_fixed <- merge(internal_fixed, degeneracy_data)
internal_fixed <- internal_fixed %>% mutate(dN=(nonsyn_count/possible_nonsyn), dS=(syn_count/possible_syn)) %>% 
  mutate(dNdS=(dN/dS))
ggplot(internal_fixed, aes(x=dN, y=dS, fill=dNdS, col=dNdS)) + 
  geom_point() + theme_bw()
labels2 <- internal_fixed %>% filter(dNdS >1)
labels2 <- labels2[!is.infinite(labels2$dNdS),]

ggplot(internal_fixed, aes(y=dN, x=dS, fill=dNdS, col=dNdS)) + 
  geom_point() + theme_bw() + 
  geom_label(data=labels2, aes(label=gene), 
             col="black",fill=alpha("white",0.2))
internal_fixed_2 <- internal_fixed %>% filter(dS > 0)
quantile(internal_fixed_2$dNdS, na.rm=T, c(0.9,0.99,0.999))
pc99_fixed <- internal_fixed_2 %>% filter(dNdS > 1.2474480)
tmp <- meta %>% filter(gene %in% pc99_fixed$gene)

## quick summary
# genes:1: STMMW_16451_1
# 2: STMMW_19031_1
# 3: STMMW_28551_1
# 4: STMMW_33691_1
# 5: STMMW_39421_1
# 6: STMMW_12961_1
# 7: STMMW_15951_1
# 8: STMMW_16821_1
# 9: STMMW_18461_1
# 10: STMMW_04181_1
# primarily pseudogenes/hypothetical genes; there is one accessory gene which is an inner membrane protein. 
# cgMLST genes: chemotaxis protein CheW; SpaM (also called invI); septum formation inhibitor Maf; 
# ubiquinone biosynthesis protein UbiJ
# nothing changed taking out those dodgy moderns

#### 2. annotate nodes with serovar label and host specific or non-specific ####
# note - if run on a slightly different dataset or tree, these node IDs will need to be updated
# host-specific
# make vectors of nodes and offspring nodes
# Typhi = N4
typhi <- c("N4")
for ( i in offspring(tr,nodeid(tr,"N4"))) {
  typhi <- rbind(typhi, nodelab(tr, i))
}
typhi
# PTA N12
paratyphiA <- c("N12")
for ( i in offspring(tr,nodeid(tr,"N12"))) {
  paratyphiA <- rbind(paratyphiA, nodelab(tr, i))
}
# both paraA and Typhi = N3

# Gallinarum: N40
gallinarum <- c("N40")
for ( i in offspring(tr,nodeid(tr,"N40"))) {
  gallinarum <- rbind(gallinarum, nodelab(tr, i))
}
#abortusovis N62
abortusovis <- c("N62")
for ( i in offspring(tr,nodeid(tr,"N62"))) {
  abortusovis <- rbind(abortusovis, nodelab(tr, i))
}
# Fulica N78
fulica <- c("N78")
for ( i in offspring(tr,nodeid(tr,"N78"))) {
  fulica <- rbind(fulica, nodelab(tr, i))
}
# abortusequi N86
abortusequi <- c("N86")
for ( i in offspring(tr,nodeid(tr,"N86"))) {
  abortusequi <- rbind(abortusequi, nodelab(tr, i))
}
# choleraesuis (modern) N129
choleraesuis_mod <- c("N129")
for ( i in offspring(tr,nodeid(tr,"N129"))) {
  choleraesuis_mod <- rbind(choleraesuis_mod, nodelab(tr, i))
}
# choleraesuis with TBS012 on N124
choleraesuis_anc <- c("N124")
for ( i in offspring(tr,nodeid(tr,"N124"))) {
  choleraesuis_anc <- rbind(choleraesuis_anc, nodelab(tr, i))
} # need to setdiff with chol modern
choleraesuis_anc <- setdiff(choleraesuis_anc,choleraesuis_mod)
#Typhisuis N140
typhisuis <- c("N140")
for ( i in offspring(tr,nodeid(tr,"N140"))) {
  typhisuis <- rbind(typhisuis, nodelab(tr, i))
}
#modern paratyphi C N157; with Tepos N156 # with TGP N149
# with KSZ N150
paratyphiC_mod <- c("N157")
for ( i in offspring(tr,nodeid(tr,"N157"))) {
  paratyphiC_mod <- rbind(paratyphiC_mod, nodelab(tr, i))
}
paratyphiC_with_tepos <- c("N156")
for ( i in offspring(tr,nodeid(tr,"N156"))) {
  paratyphiC_with_tepos <- rbind(paratyphiC_with_tepos, nodelab(tr, i))
}
paratyphiC_with_ksz <- c("N150")
for ( i in offspring(tr,nodeid(tr,"N150"))) {
  paratyphiC_with_ksz <- rbind(paratyphiC_with_ksz, nodelab(tr, i))
}
paratyphiC_with_tgp <- c("N149")
for ( i in offspring(tr,nodeid(tr,"N149"))) {
  paratyphiC_with_tgp <- rbind(paratyphiC_with_tgp, nodelab(tr, i))
}
# do setdiff at the end
paratyphiC_with_tgp <- setdiff(paratyphiC_with_tgp,paratyphiC_with_ksz)
paratyphiC_with_ksz <- setdiff(paratyphiC_with_ksz,paratyphiC_with_tepos)
paratyphiC_with_tepos <- setdiff(paratyphiC_with_tepos,paratyphiC_mod)
# host generalist
#Infantis N23
infantis <- c("N22")
for ( i in offspring(tr,nodeid(tr,"N22"))) {
  infantis <- rbind(infantis, nodelab(tr, i))
}
# Tallahassee N52
tallahassee <- c("N52")
for ( i in offspring(tr,nodeid(tr,"N52"))) {
  tallahassee <- rbind(tallahassee, nodelab(tr, i))
}
# Bispebjerg N69
bispebjerg <- c("N69")
for ( i in offspring(tr,nodeid(tr,"N69"))) {
  bispebjerg <- rbind(bispebjerg, nodelab(tr, i))
}
# Berta N29
berta <- c("N29")
for ( i in offspring(tr,nodeid(tr,"N29"))) {
  berta <- rbind(berta, nodelab(tr, i))
}
# Birkenhead N97
birkenhead <- c("N97")
for ( i in offspring(tr,nodeid(tr,"N97"))) {
  birkenhead <- rbind(birkenhead, nodelab(tr, i))
}
## ancients
# need to setdiff with other vectors
# N41: birkenhead with ancients 
bkhd_ancs <- c("N97")
for ( i in offspring(tr,nodeid(tr,"N97"))) {
  bkhd_ancs <- rbind(bkhd_ancs, nodelab(tr, i))
}
# N107 = bronze age lineage 
ba_ancs <- c("N107")
for ( i in offspring(tr,nodeid(tr,"N107"))) {
  ba_ancs <- rbind(ba_ancs, nodelab(tr, i))
}
# N120 KNC clade with XBQM and BPN
knc_cl <- c("N120")
for ( i in offspring(tr,nodeid(tr,"N120"))) {
  knc_cl <- rbind(knc_cl, nodelab(tr, i))
}
# N51- all ancients except birkenhead, need to setdiff with prev
ancs_except_bkhd <- c("N106")
for ( i in offspring(tr,nodeid(tr,"N106"))) {
  ancs_except_bkhd <- rbind(ancs_except_bkhd, nodelab(tr, i))
}

ancs_except_bkhd <- setdiff(ancs_except_bkhd,c(abortusequi,abortusovis,berta,birkenhead,bispebjerg,choleraesuis_anc,choleraesuis_mod,
                                               fulica,gallinarum,infantis,paratyphiA,paratyphiC_mod,paratyphiC_with_ksz,paratyphiC_with_tepos,
                                               paratyphiC_with_tgp,tallahassee,typhi,typhisuis,bkhd_ancs,ba_ancs,knc_cl))

# N40 - all ancients incl bkhd need to setdiff. 
ancs_incl_bkhd <- c("N95")
for ( i in offspring(tr,nodeid(tr,"N95"))) {
  ancs_incl_bkhd <- rbind(ancs_incl_bkhd, nodelab(tr, i))
}
ancs_incl_bkhd <- setdiff(ancs_incl_bkhd,c(abortusequi,abortusovis,berta,birkenhead,bispebjerg,choleraesuis_anc,choleraesuis_mod,
                                           fulica,gallinarum,infantis,paratyphiA,paratyphiC_mod,paratyphiC_with_ksz,paratyphiC_with_tepos,
                                           paratyphiC_with_tgp,tallahassee,typhi,typhisuis,bkhd_ancs,ba_ancs,knc_cl,ancs_except_bkhd))


## for each vector, annotate the table
# don't forget to also include:
# shared paraA and Typhi = N3 (included in host-specific classification)
#shared abequi and fulica is N156
df <- df %>% mutate(serovar=if_else(node %in% abortusequi,"abortusequi",
if_else(node %in% abortusovis,"abortusovis",
if_else(node %in% choleraesuis_anc,"choleraesuis_anc",
if_else(node %in% choleraesuis_mod,"choleraesuis_mod",
if_else(node %in% fulica,"fulica",
if_else(node %in% gallinarum,"gallinarum",
if_else(node %in% paratyphiA,"paratyphiA",
if_else(node %in% paratyphiC_mod,"paratyphiC_mod",
if_else(node %in% paratyphiC_with_ksz,"paratyphiC_with_ksz",
if_else(node %in% paratyphiC_with_tepos,"paratyphiC_with_tepos",
if_else(node %in% paratyphiC_with_tgp,"paratyphiC_with_tgp",
if_else(node %in% typhi,"typhi",
if_else(node %in% typhisuis,"typhisuis",
if_else(node == "N3","shared_paraA_typhi",
if_else(node=="N77","shared_abequi_fulica", 
if_else(node %in% berta,"berta",
if_else(node %in% birkenhead,"birkenhead",
if_else(node %in% bispebjerg,"bispebjerg",
if_else(node %in% infantis,"infantis",
if_else(node %in% tallahassee,"tallahassee",
if_else(node %in% bkhd_ancs,"bkhd_ancs",
if_else(node %in% ba_ancs,"ba_ancs",
if_else(node %in% knc_cl,"knc_cl",
if_else(node %in% ancs_except_bkhd,"ancs_except_bkhd",
if_else(node %in% ancs_incl_bkhd,"ancs_incl_bkhd","unassigned_ancestral"))))))))))))))))))))))))))


# for labelling as host-adapted or not
host_adapt <- c("abortusequi","abortusovis","choleraesuis_anc","choleraesuis_mod","fulica",
                "gallinarum","paratyphiA","paratyphiC_mod","paratyphiC_with_ksz","paratyphiC_with_tepos",
                "paratyphiC_with_tgp","typhi","typhisuis","shared_paraA_typhi","shared_abequi_fulica")
# 
host_gen <- c("berta","birkenhead","bispebjerg","infantis","tallahassee","unassigned_ancestral")
# 
older_ancients <- c("bkhd_ancs","ba_ancs","knc_cl","ancs_except_bkhd","ancs_incl_bkhd")
df <- df %>% mutate(adapted_generalist=if_else(serovar %in% host_adapt,"host_adapted",
                                               if_else(serovar %in% older_ancients,"ancients_unassigned",
                                                       "host_generalist")))

# save df 
fwrite(df, "snppar_annotations.with_fixed_serovar_and_adapted_annotations.csv")
#### end annotation ####
#### 3. plots for node-level summaries ####
df
# subset to internal
df_internal <- df %>% filter(is_internal==TRUE)
## check revertant histogram
hist(df_internal$total_rev_nodes/df_internal$n_descendant_nodes)
df_internal_max5pc_rev <- df_internal %>% filter(total_rev_nodes/n_descendant_nodes < 0.05)

df_internal_max5pc_rev <- df_internal_max5pc_rev %>% mutate(nonsyn=if_else(impact %in% c("high","missense"),1,0)) %>% 
  mutate(syn=if_else(impact %in% c("synonymous"),1,0))
node_summary_internal <- df_internal_max5pc_rev %>% group_by(node) %>% mutate(nonsyn_count=sum(nonsyn), syn_count=sum(syn)) %>% 
  select(node, serovar, adapted_generalist, nonsyn_count, syn_count) %>% unique()
node_summary_internal <- node_summary_internal %>% mutate(dnds=nonsyn_count/syn_count)
# plot dN/dS by node type
p1_violin <- ggplot(node_summary_internal, aes(y=dnds, x=adapted_generalist)) + geom_violin() + geom_point() + theme_bw() +
  ggrepel::geom_label_repel(aes(label=serovar)) + ggtitle("max 5% revertant internal mutations by node label") + ylim(0,10)

p2_violin <- node_summary_internal %>% filter(syn_count > 1) %>% ggplot(aes(y=dnds, x=adapted_generalist)) + geom_violin() + geom_point() + theme_bw() +
  ggrepel::geom_label_repel(aes(label=serovar)) + ggtitle("max 5% revertant internal mutations by node label; minimum 2 synonymous mutations at a node")

p1_box <- ggplot(node_summary_internal, aes(y=dnds, x=adapted_generalist)) + geom_boxplot() + geom_point() + theme_bw() +
  ggrepel::geom_label_repel(aes(label=serovar)) + ggtitle("max 5% revertant internal mutations by node label") + ylim(0,10)

p2_box <- node_summary_internal %>% filter(syn_count > 1) %>% ggplot(aes(y=dnds, x=adapted_generalist)) + geom_boxplot() + geom_point() + theme_bw() +
  ggrepel::geom_label_repel(aes(label=serovar)) + ggtitle("max 5% revertant internal mutations by node label; minimum 2 synonymous mutations at a node")

#### 4. plot with tree data ####
# function to add labels to tree
add_labels <- function(tree, data) {
  missing_numbers <- c()
  
  for (i in 1:nrow(data)) {
    label_name <- data$node[i]
    label_number <- data$dnds[i]
    print(label_name)
    print(label_number)
    if (label_name %in% tree$tip.label) {
      tip_index <- match(label_name, tree$tip.label)
      #tree$tip.label[tip_index] <- paste(label_name, label_number, sep = ":")
      tree$dnds[tip_index] <- label_number
    } 
    else if (label_name %in% tree$node.label) {
      node_index <- match(label_name, tree$node.label)
      #tree$node.label[node_index] <- paste(label_name, label_number, sep = ":")
      tree$dnds[length(tree$tip.label) + node_index] <- label_number
    } 
    else {
      missing_numbers <- c(missing_numbers, label_name)
    }
  }
  
  return(list(tree = tree, missing_numbers = missing_numbers))
}
# node_summary_internal and tr - what you need here 
# subset to min 2 syn muts?
min2_syn_internal <- node_summary_internal %>% filter(syn_count > 1)
# create a dataframe with tip labels also
all_nodelabs <- tr %>% as.tibble() %>% select(label) %>% unique()
diff_nodes <- setdiff(all_nodelabs$label, min2_syn_internal$node)
dummy_data <- min2_syn_internal %>% select(node, dnds)
dummy_data$dnds <- as.character(dummy_data$dnds) # temporary
tmp <- c("node","dnds")
names(tmp) <- c("node","dnds")
for (i in 1:length(diff_nodes)) {
  tmp[1] <- diff_nodes[i]
  tmp[2] <- "NA"
  #  print(tmp)
  dummy_data <- rbind(dummy_data,tmp)
}
dummy_data$dnds <- as.numeric(dummy_data$dnds) 
dummy_data
nrow(dummy_data)
# Modify the tree object to include labels and numbers
# test
result <- add_labels(tr, dummy_data)
tr_new <- result$tree
tr_new$dnds <- as.numeric(tr_new$dnds)
tr_new$dnds
# use a log10 transformation so 1 = 0
coolwarm_hcl <- colorspace::diverging_hcl(11,
                                          h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
ggtree(tr_new, layout="circular", branch.length="none") + 
  geom_nodepoint(aes(subset=!isTip, fill=(log10(tr_new$dnds)), col=(log10(tr_new$dnds)))) + 
  scale_fill_gradient2(low = coolwarm_hcl[1], high =  coolwarm_hcl[11], mid = "white", na.value = "darkgrey")+ 
  scale_color_gradient2(low = coolwarm_hcl[1], high =  coolwarm_hcl[11],mid = "white",na.value="darkgrey")+ # remove na points later
  geom_cladelab(node = nodeid(tr,"N129"), "Modern Choleraesuis")+
  geom_cladelab(node = nodeid(tr,"N157"), "Modern Paratyphi C")+
  geom_cladelab(node = nodeid(tr,"N140"), "Modern Typhisuis")+
  geom_cladelab(node = nodeid(tr,"N4"), "Typhi")+
  geom_cladelab(node = nodeid(tr,"N107"), "BA lineage")+
  geom_cladelab(node = nodeid(tr,"N97"), "Modern Birkenhead")+
  geom_cladelab(node = nodeid(tr,"N12"), "Paratyphi A")+
  geom_cladelab(node = nodeid(tr,"N22"), "Infantis")+
  geom_cladelab(node = nodeid(tr,"N40"), "Gallinarum")+
  geom_cladelab(node = nodeid(tr,"N52"), "Tallahassee")+ 
  geom_cladelab(node = nodeid(tr,"N62"), "Abortusovis")+
  geom_cladelab(node = nodeid(tr,"N69"), "Bispebjerg")+
  geom_cladelab(node = nodeid(tr,"N78"), "Fulica-like")+
  geom_cladelab(node = nodeid(tr,"N86"), "Abortusequi")+
  geom_cladelab(node = nodeid(tr,"N29"), "Berta")
# plot without a log10 transformation
ggtree(tr_new, layout="circular", branch.length="none") + 
  geom_nodepoint(aes(subset=!isTip, fill=(tr_new$dnds), col=(tr_new$dnds))) + 
  scale_fill_gradient2(low = coolwarm_hcl[1], high =  coolwarm_hcl[11], mid = "white", na.value = "darkgrey")+ 
  scale_color_gradient2(low = coolwarm_hcl[1], high =  coolwarm_hcl[11],mid = "white",na.value="darkgrey")+ # remove na points later
  geom_cladelab(node = nodeid(tr,"N129"), "Modern Choleraesuis")+
  geom_cladelab(node = nodeid(tr,"N157"), "Modern Paratyphi C")+
  geom_cladelab(node = nodeid(tr,"N140"), "Modern Typhisuis")+
  geom_cladelab(node = nodeid(tr,"N4"), "Typhi")+
  geom_cladelab(node = nodeid(tr,"N107"), "BA lineage")+
  geom_cladelab(node = nodeid(tr,"N97"), "Modern Birkenhead")+
  geom_cladelab(node = nodeid(tr,"N12"), "Paratyphi A")+
  geom_cladelab(node = nodeid(tr,"N22"), "Infantis")+
  geom_cladelab(node = nodeid(tr,"N40"), "Gallinarum")+
  geom_cladelab(node = nodeid(tr,"N52"), "Tallahassee")+ 
  geom_cladelab(node = nodeid(tr,"N62"), "Abortusovis")+
  geom_cladelab(node = nodeid(tr,"N69"), "Bispebjerg")+
  geom_cladelab(node = nodeid(tr,"N78"), "Fulica-like")+
  geom_cladelab(node = nodeid(tr,"N86"), "Abortusequi")+
  geom_cladelab(node = nodeid(tr,"N29"), "Berta")
# min 5 doesn't make this clearer 
#### end tree plots ####
#### 5. Convergent mutations plotting ####
df_internal_max5pc_rev <- df_internal_max5pc_rev %>% mutate(varid=paste0(gene,"_", pos, "_", impact))%>% 
  mutate(varid2=paste0(gene,"_", impact))

## start by plotting overlaps between different host-specialist and host-generalist serovars 

# all non-synonymous at gene level
ns_gene_level <- list(Generalist=unique(df_internal_max5pc_rev[adapted_generalist=="host_generalist" & nonsyn==1 ,varid2]),
                      ParatyphiC=unique(df_internal_max5pc_rev[serovar=="paratyphiC_mod" & nonsyn==1 ,varid2]),
                      Typhi=unique(df_internal_max5pc_rev[serovar=="typhi" & nonsyn==1 ,varid2]),
                      Abortusequi =unique(df_internal_max5pc_rev[serovar=="abortusequi" & nonsyn==1 ,varid2]),
                      Abortusovis =unique(df_internal_max5pc_rev[serovar=="abortusovis" & nonsyn==1 ,varid2]),
                      Fulica = unique(df_internal_max5pc_rev[serovar=="fulica" & nonsyn==1 ,varid2]),
                      Gallinarum =unique(df_internal_max5pc_rev[serovar=="gallinarum" & nonsyn==1 ,varid2]),
                      ParatyphiA=unique(df_internal_max5pc_rev[serovar=="paratyphiA" & nonsyn==1 ,varid2]), 
                      BA=unique(df_internal_max5pc_rev[serovar=="ba_ancs" & nonsyn==1 ,varid2]),
                      ParatyphiC_with_ancs=unique(df_internal_max5pc_rev[serovar %in% c("paratyphiC_with_ksz","paratyphiC_with_tepos","paratyphiC_with_tgp") & nonsyn==1 ,varid2]),
                      Choleraesus_w_ancs_and_mod=unique(df_internal_max5pc_rev[serovar %in% c("choleraesuis_mod","choleraesuis_anc") & nonsyn==1 ,varid2]),
                      Typhisuis=unique(df_internal_max5pc_rev[serovar=="typhisuis" & nonsyn==1 ,varid2]))
## this is a mess
plot_upset(Venn(ns_gene_level), 
           nintersects = 77)
# all synonymous at gene level # incomprehensible
# all pseudo at gene level
pseudo_gene_level <- list(Generalist=unique(df_internal_max5pc_rev[adapted_generalist=="host_generalist" & impact=="high" ,varid2]),
                          ParatyphiC=unique(df_internal_max5pc_rev[serovar=="paratyphiC_mod" & impact=="high" ,varid2]),
                          Typhi=unique(df_internal_max5pc_rev[serovar=="typhi" & impact=="high" ,varid2]),
                          Abortusequi =unique(df_internal_max5pc_rev[serovar=="abortusequi" & impact=="high" ,varid2]),
                          Abortusovis =unique(df_internal_max5pc_rev[serovar=="abortusovis" & impact=="high" ,varid2]),
                          Fulica = unique(df_internal_max5pc_rev[serovar=="fulica" & impact=="high" ,varid2]),
                          Gallinarum =unique(df_internal_max5pc_rev[serovar=="gallinarum" & impact=="high" ,varid2]),
                          ParatyphiA=unique(df_internal_max5pc_rev[serovar=="paratyphiA" & impact=="high" ,varid2]), 
                          BA=unique(df_internal_max5pc_rev[serovar=="ba_ancs" & impact=="high" ,varid2]),
                          ParatyphiC_with_ancs=unique(df_internal_max5pc_rev[serovar %in% c("paratyphiC_with_ksz","paratyphiC_with_tepos","paratyphiC_with_tgp") &impact=="high" ,varid2]),
                          Choleraesus_w_ancs_and_mod=unique(df_internal_max5pc_rev[serovar %in% c("choleraesuis_mod","choleraesuis_anc") & impact=="high" ,varid2]),
                          Typhisuis=unique(df_internal_max5pc_rev[serovar=="typhisuis" & nonsyn==1 ,varid2]))

plot_upset(Venn(pseudo_gene_level), 
           nintersects = 77)

# here, you don't have the common ancestor of typhi + paratyphi A assignments in (or abortusequs and fulica), 
# which then I think cuts out a chunk of data that you'd actually be interested in
# non-synonymous overlaps at position level
ns_var_level <- list(Generalist=unique(df_internal_max5pc_rev[adapted_generalist=="host_generalist" & nonsyn==1 ,varid]),
                     ParatyphiC=unique(df_internal_max5pc_rev[serovar=="paratyphiC_mod" & nonsyn==1 ,varid]),
                     Typhi=unique(df_internal_max5pc_rev[serovar=="typhi" & nonsyn==1 ,varid2]),
                     Abortusequi =unique(df_internal_max5pc_rev[serovar=="abortusequi" & nonsyn==1 ,varid]),
                     Abortusovis =unique(df_internal_max5pc_rev[serovar=="abortusovis" & nonsyn==1 ,varid]),
                     Fulica = unique(df_internal_max5pc_rev[serovar=="fulica" & nonsyn==1 ,varid]),
                     Gallinarum =unique(df_internal_max5pc_rev[serovar=="gallinarum" & nonsyn==1 ,varid]),
                     ParatyphiA=unique(df_internal_max5pc_rev[serovar=="paratyphiA" & nonsyn==1 ,varid]), 
                     BA=unique(df_internal_max5pc_rev[serovar=="ba_ancs" & nonsyn==1 ,varid]),
                     ParatyphiC_with_ancs=unique(df_internal_max5pc_rev[serovar %in% c("paratyphiC_with_ksz","paratyphiC_with_tepos","paratyphiC_with_tgp") & nonsyn==1 ,varid]),
                     Choleraesus_w_ancs_and_mod=unique(df_internal_max5pc_rev[serovar %in% c("choleraesuis_mod","choleraesuis_anc") & nonsyn==1 ,varid]),
                     Typhisuis=unique(df_internal_max5pc_rev[serovar=="typhisuis" & nonsyn==1 ,varid]))
plot_upset(Venn(ns_var_level), 
           nintersects = 72)
#### end ####

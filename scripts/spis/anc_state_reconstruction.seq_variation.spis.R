rm(list=ls())
# library imports
library(phytools)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
library(treeio)
library(ggtree)
library(cowplot)
#### Anc state reconstruction #### 
# 1. ancestral state reconstruction with phytools
# need tree + phenotype table
# first - midpoint root tree
ml_tree_paraC <- 'scaffold_paraC_mq37_mincov3_excl_gff.geno05.iqtree.treefile'
tr <- ape::read.tree()
tr<- midpoint_root(tr)
# save tree
#write.tree(tr,"scaffold_paraC_mq37_mincov3_excl_gff.geno05.iqtree.midp_root.treefile")
paraC_tr <- read.newick("scaffold_paraC_mq37_mincov3_excl_gff.geno05.iqtree.midp_root.treefile",
                        node.label='support')
pres_absence <- 'SPIs_pres_absence.tab' # 1 if "present" based on breadth of coverage, 0 if absent
spis <- read.csv(pres_absence,row.names = 1, sep = "\t")
spis
# split to individual islands for asr
sp1 <- as.matrix(spis)[,1]
sp2 <- as.matrix(spis)[,2]
sp3 <- as.matrix(spis)[,3]
sp4 <- as.matrix(spis)[,4]
sp5 <- as.matrix(spis)[,5]
sp6 <- as.matrix(spis)[,6]
sp7 <- as.matrix(spis)[,7]
sp8 <- as.matrix(spis)[,8]
sp9 <- as.matrix(spis)[,9]
sp10 <- as.matrix(spis)[,10]

#SPI-7 is only one which is variable - there is some missingness in modern birkenhead, but not so much
is.rooted(paraC_tr) # sanity check

ggtree(paraC_tr) + geom_tiplab() + geom_nodelab(aes(label=node))

v <- as.vector(sp7)
names(v) <- names(sp7)

cols<-setNames(palette()[1:length(unique(sp7))],sort(unique(sp7)))
# ancestral state reconstruction with diff models
# no bootstrapping 
fitER<-ace(sp7,paraC_tr@phylo,model="ER",type="discrete")
fitARD <- ace(sp7,paraC_tr@phylo,model="ARD",type="discrete")
fitER
fitARD
cbind(fitER$lik.anc,fitARD$lik.anc) # double check these are basically the same

ancstate <- fitARD$lik.anc
#names <-  rownames(ancstate)
ancstate <- data.frame(ancstate)
ancstate$names <- rownames(ancstate)
colnames(ancstate) <- c('absent','present','node')
p <- ggtree(paraC_tr) + geom_tiplab() 
pies <- nodepie(ancstate, cols = 1:2)
cols2 <- cols
names(cols2) <- c('absent','present')
pies <- lapply(pies, function(g) g+scale_fill_manual(values = cols2))
ard_tr <- p + geom_inset(pies, width = .025, height = .025) 

cowplot::plot_grid(er_tr,ard_tr,ncol=1,labels=c("ER","ARD"))

ggtree(paraC_tr,branch.length="none") + geom_tiplab() + geom_inset(pies, width = .025, height = .025)


## simulate single stochastic character map using empirical Bayes method
mtree<-make.simmap(paraC_tr@phylo,sp7,model="ARD",nsim=500) # use 500 sims
# save this, i don't want to have to re-run it every time
saveRDS(mtree,"spi7_anc_reconstruction.make_simmap.500_sims.ARD_model.rds")
mtree <- readRDS("spi7_anc_reconstruction.make_simmap.500_sims.ARD_model.rds")

names(cols2) <- c(0,1) # make sure colours saved properly
# plot quickly
par(mfrow=c(10,10))
null<-sapply(mtree,plot,colors=cols2,lwd=1,ftype="off")

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(paraC_tr@phylo)),fsize=0.8)

pd<-summary(mtree,plot=FALSE)
pd
plot(pd,ftype="i")


## Plot sim results properly 
node_proportions <- data.frame(describe.simmap(mtree)$ace)

node_proportions$names <- rownames(node_proportions)
names(node_proportions)

colnames(node_proportions) <- c('absent','present','node')
# remove any tips from this dataframe -> have a seperate data frame for tip and nodes
nps <- node_proportions %>% filter(!node %in% paraC_tr@phylo$tip.label)
tps <- node_proportions %>% filter(node %in% paraC_tr@phylo$tip.label)

p <- ggtree(paraC_tr) + geom_tiplab() 
p
pies <- nodepie(nps, cols = 1:2)
names(cols2) <- c('absent','present')
pies2 <- lapply(pies, function(g) g + scale_fill_manual(values = cols2))

p + geom_inset(pies2, width = .025, height = .025) 
#### GRY004 ####
# 2. Heterozygosity - what is causing GRY004 result
library(cowplot)
# read in full het results
het_gry <- '/path/to/qc/individual_csvs/GRY004.csv.gz'
het <- fread(het_gry)
# filter for sites where you actually see heterozygosity
het <- het %>% filter(Depth > 0) %>% filter(Minor_Support > 0)

ggplot(het, aes(x=Minor_Support)) + geom_histogram() + 
  facet_wrap(~QualFilt+MutType, scales="free_y") +
  theme_cowplot()

het_qualfilt <- het %>% filter(QualFilt != "LowQual")
ggplot(het_qualfilt, aes(x=Minor_Support)) + geom_histogram() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot()
ggplot(het_qualfilt, aes(x=Minor_Support)) + geom_density() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot()
summary(het_qualfilt$Minor_Support)
ggplot(het_qualfilt, aes(x=Depth)) + geom_density() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot()
ggplot(het_qualfilt, aes(x=Minor_Support,y=Depth)) + geom_density_2d() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot()
ggplot(het_qualfilt, aes(x=Minor_Support,y=Depth)) + geom_density_2d_filled() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot()
# for transitions, there are two "peaks" - at a little less than 0.1 and at around 0.2
# 0.2 tend to be higher depth overall as well 
# for transversions- this is only at the 0,1 mark, with much lower density around the 0.2 mark

# distribution throughout genome
ggplot(het_qualfilt, aes(x=Position,y=Minor_Support*Depth)) + geom_point() + theme_cowplot()+ 
  facet_wrap(~MutType, scales="free_y") 
# not a clear pattern. 

# how does het compare to coverage reductions at genomic islands? 
hetcov <- 'path/to/qualimap_data_combined/coverage_across_reference_spi7_positive.no_sd.tab'
cov <- fread(hetcov)

head(cov)
colnames(cov) <- c("Position","Depth","ID")
#order by phylo position
cov$ID <- factor(cov$ID, levels=c("GRY004","GRY001","KPI001","KKN099","AKY001",
                                  "KNC041","KNC056","KNC095",
                                  "KSZ005","KNL026","HDL004","MGL008",
                                  "TAV007","HGH1600","Tepos_10",
                                  "Tepos_14","Tepos_35","Tepos_37"))

ggplot(cov, aes(x=Position,y=Depth)) + geom_line() + 
  geom_vline(xintercept=4409511, linetype="dotted",col="red") +
  geom_vline(xintercept=4543148, linetype="dotted",col="red") +
  facet_wrap(~ID, scales="free_y", ncol = 2) + 
  cowplot::theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )

# for each individual - mean coverage at SPI vs mean coverage overall
# label windows as SPI-7 vs non spi 7
cov2 <- cov %>% mutate(spi=if_else((Position > 4409511 & Position < 4543148),'spi7','no')) 

cov2 <- cov2 %>% group_by(ID,spi) %>% summarise(mean_dp=mean(Depth),sd_dp=sd(Depth),median_dp=median(Depth))
tmp1 <- cov2 %>% filter(spi=='spi7')
tmp2 <- cov2 %>% filter(spi=='no')
colnames(tmp1) <- c('ID','spi','mean_dp_spi7','sd_dp_spi7','median_dp_spi7')
colnames(tmp2) <- c('ID','spi','mean_dp_gx','sd_dp_gx','median_dp_gx')

tmp1$spi  <- NULL
tmp2$spi <- NULL
cov2 <- merge(tmp1,tmp2, by="ID")
cov2 <- cov2 %>% mutate(mean_ratio=mean_dp_spi7/mean_dp_gx,
               median_ratio=median_dp_spi7/median_dp_gx)
library(cowplot)
ggplot(cov2, aes(x=mean_ratio)) + geom_histogram() + theme_cowplot()
ggplot(cov2, aes(x=median_ratio))+ geom_histogram() + theme_cowplot()

p1 <- ggplot(cov, aes(x=Position,y=Depth)) + geom_line() + geom_point() +
  xlim(4409511,4543148) + 
  facet_wrap(~ID, scales="free_y", ncol = 1) + 
  geom_hline(data=cov2, aes( yintercept=mean_dp_gx),linetype="dotted",colour="red") +
  cowplot::theme_cowplot()

plt2 <- ggplot(cov, aes(x=Position,y=Depth)) + geom_line() + geom_point() +
  xlim(4409511,4543148) + 
  facet_wrap(~ID, scales="free_y", ncol = 1) + 
  geom_hline(data=cov2, aes( yintercept=mean_dp_gx),linetype="dotted",colour="red") +
  cowplot::theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
  )

plt2


# restricting to SPI7 region, what does minor support density plot look like for GRY genomes

het_qualfilt_spi7 <- het_qualfilt %>% filter(Position > 4409511)%>% filter(Position < 4543148)
# nothing

full_genome_density <- ggplot(het_qualfilt, aes(x=Minor_Support,y=Depth)) + geom_density_2d_filled() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot() + guides(fill="none")

cov_ratios <- ggplot(cov2, aes(x=mean_ratio)) + geom_histogram() + theme_cowplot() + 
  xlab("Mean Depth SPI-7/Mean Depth Genome") + ylab("# Genomes") +
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5))

plot_grid(cov_ratios,full_genome_density, ncol = 1, labels=c("A","B"))
ggplot(cov2, aes(x=median_ratio))+ geom_histogram() + theme_cowplot()
nrow(het_qualfilt)
ggplot(het_qualfilt_spi72, aes(x=Minor_Support,y=Depth)) + geom_density_2d_filled() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot() + ylim(0,10)

# need to exclude problematic positions I think 
# position in paratyphi C reference: 4503961..4519253
to_excl <- fread('data/paraC_mask.txt') # positions from GFF which should be excluded 
het_qualfilt <- het_qualfilt%>% filter(!Position %in% to_excl$V1)
het_qualfilt_spi7 <- het_qualfilt %>% filter(Position > 4503961)%>% filter(Position < 4519253)
ggplot(het_qualfilt, aes(x=Minor_Support,y=Depth)) + geom_density_2d_filled() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot() + guides(fill="none")

ggplot(het_qualfilt_spi7, aes(x=Minor_Support,y=Depth)) + geom_density_2d_filled() + 
  facet_wrap(~MutType, scales="free_y") +
  theme_cowplot() + guides(fill="none")


#### Plotting coverage vs gene types #### 

# gene type info 
spi6_df <- fread('data/SPI6_features_to_plot.tsv',sep='\t')
spi7_df <- fread('data/SPI7_features_to_plot.tsv',sep='\t')

# coverage dfs 
spi_genecov <- fread('/path/to/genecov/table')
colnames(spi_genecov) <- c('chromosome','start','end','id_bedtools','nreads','nbases_covered','length','frac_covered','sample_id')
spi7_genecov <- spi_genecov[grepl("SPI-7",spi_genecov$id_bedtools)]
spi6_genecov <- fread('/path/to/genecov/table')
colnames(spi6_genecov) <- c('chromosome','start','end','id_bedtools','nreads','nbases_covered','length','frac_covered','sample_id')

head(spi6_genecov)
spi6_genecov_with_meta <- merge(spi6_genecov, spi6_df, by="end")
# colour=classification),
# linewidth=3) + theme_cowplot()

#order IDs by position in main figure 
id_ord_ptc <- c("Schwarzengrund-CVM19633","BareillyCFSAN000189","Bovismorbificans3114",
                "EnteritidisP125109","GRY004",
                'SME005','SJL020','SJL019','IV3002','IKI003','BOY008',
                'S1L002','HGC040','HGC004','KKN099','AKY001','KPI001','GRY001','H1I004',
                'MAJ022','SUA004','KO1037','SAP004','HO1003','NEP019','NEP021','NEP018',
                'CHL003','XBQM20','KNC181','KNC056','KNC095','KNC041','BPN005','XBQM90',
                'CPA002','SIA006','TBS012','ETR001','BRC034','BRC017',
                'MGL011','CholeraesuisSC-B67','KNC139',
                'TGP001','KSZ005','KNL026','MGL008','HDL004',
                'GTH003','TAV007','KNL019','Ragna','HGH1600',
                'Tepos_10','Tepos_14','Tepos_35',
                'Tepos_37','ParatyphiC-RKS4594','TyphiTy2')

spi6_genecov$sample_id <- factor(spi6_genecov$sample_id, levels=id_ord_ptc)
p1 <- ggplot(spi6_genecov) + geom_segment(aes(x=start,xend=end,colour = frac_covered,y=sample_id), linewidth=3) +
  scale_color_viridis_c()+ theme_cowplot()
# add a facet for gene type
p2 <- ggplot(spi6_df) + geom_segment(aes(x=start,xend=end, colour = `gene type`, y=strand),linewidth=4) + theme_cowplot()
plot_grid(p1,p2, ncol=1, rel_heights=c(0.9,0.2))
p3 <-  ggplot(spi6_df) + geom_segment(aes(x=start,xend=end, colour = `brite hierarchies`, y=strand),linewidth=4) + theme_cowplot()
plot_grid(p1,p3, ncol=1, rel_heights=c(0.9,0.2))
p3
unique(spi6_df$`brite hierarchies`)
spi6_df<- spi6_df %>% mutate(class_summary=if_else(`brite hierarchies` %in% c("annotated as hypothetical protein kegg","not annotated kegg","conserved hypothetical protein"),"hypothetical",
                      if_else(`brite hierarchies` %in% c("02044 Secretion system [BR:sty02044]" ,"02044 Secretion system [BR:sty02044]; 02035 Bacterial motility proteins [BR:sty02035]"),"secretion_sys",
                              if_else(`brite hierarchies` %in% c("rhs_family_protein","02048 Prokaryotic defense system [BR:sty02048]"),"prokaryotic_defence",
                                      if_else(`brite hierarchies` %in% c("probable secreted protein","safA lipoprotein","outer membrane adhesin","putative fimbrial protein"),"external/surface",
                                              if_else(`brite hierarchies` %in% c("fimbrial chaperone","fimbrial structural subunit", "putative fimbrial subunit","outer membrane fimbrial usher protein"),"fimbrial proteins",
                                                      if_else(`brite hierarchies`%in% c("annotated as pseudogene in KEGG","putative membrane component kegg"),"pseudo",`brite hierarchies`)))))))# %>%


p4 <- ggplot(spi6_df) + geom_segment(aes(x=start,xend=end, colour = class_summary, y=strand),linewidth=4) + theme_cowplot()
p4
plot_grid(p1,p4, ncol=1, rel_heights=c(0.9,0.2))

# also do something for spi7
spi7_genecov<- spi7_genecov %>% filter(sample_id %in% id_ord_ptc)
spi7_genecov$sample_id <- factor(spi7_genecov$sample_id, levels=id_ord_ptc)
q1 <- ggplot(spi7_genecov) + geom_segment(aes(x=start,xend=end,colour = frac_covered,y=sample_id), linewidth=3) +
  scale_color_viridis_c()+ theme_cowplot()
q1
spi7_df
q2 <- ggplot(spi7_df) + 
  geom_segment(aes(x=start,xend=end, colour = classification, y=strand),linewidth=4) + theme_cowplot()
plot_grid(q1,q2, ncol=1, rel_heights=c(0.9,0.2))

#### Sequence variation from pangenome alignment ####

# genotype tables from VCFs (generated using vcf-to-tab)
tab1 <- 'SPI6_regions.all_ancients.moderns_subset.ac1.tab'
tab2 <- 'SPI7_regions.all_ancients.moderns_subset.ac1.tab'
gt_tab_spi6 <- fread(tab1) 
gt_tab_spi7 <- fread(tab2)

# get table of genotypes for loci of interest.

#SPI-7: pil and Vi locus

pil_locus_tags <- c('T_RS21585_1','T_RS21590_1','STY4543_1','T_RS21610_1',
                    'T_RS21615_1','T_RS21620_1','STY4547_1','SPAB_RS21825_1',
                    'T_RS21635_1')

vi_locus_tags <- c('T_RS22115_1','T_RS22120_1','T_RS22125_1','T_RS22130_1',
                   'T_RS22135_1','T_RS22140_1','STY4659_1','T_RS22150_1',
                   'T_RS22155_1', 'T_RS22160_1')

# SPI-6 - saf and tcf locus tags
# note that first three saf locus tags are all safA
saf_locus_tags <- c('T_RS13015_1','SNSL254_RS02590_1','STMMW_03091_1','STMMW_03101_1',
                    'STMMW_03111_1','STMMW_03121_1')

tcf_locus_tags <- c('T_RS12955_1','T_RS12950_1','T_RS12945_1','T_RS12940_1')

saf_gts <- gt_tab_spi6 %>% filter(`#CHROM` %in% saf_locus_tags)
tcf_gts <- gt_tab_spi6 %>% filter(`#CHROM` %in% tcf_locus_tags)
tcf_gts
colnames(saf_gts)
# get IDs to keep = those in paraC tree (only)
tr <- ape::read.tree('scaffold_paraC_mq37_mincov3_excl_gff.geno05.iqtree.midp_root.treefile')

ids_to_keep <- tr$tip.label
saf_gt_tab <- saf_gts %>% select(all_of(c("#CHROM","POS","REF",ids_to_keep)))
tcf_gt_tab <- tcf_gts %>% select(all_of(c("#CHROM","POS","REF",ids_to_keep)))
vi_gts <- gt_tab_spi7 %>% filter(`#CHROM` %in% vi_locus_tags)
pil_gts <- gt_tab_spi7 %>% filter(`#CHROM` %in% pil_locus_tags)

vi_gt_tab <- vi_gts %>% select(all_of(c("#CHROM","POS","REF",ids_to_keep)))
pil_gt_tab <- pil_gts %>% select(all_of(c("#CHROM","POS","REF",ids_to_keep)))

# add a column for gene ID 
unique(vi_gts$`#CHROM`)
tsl_vi <- data.frame(cbind(c('T_RS22115_1','T_RS22120_1','T_RS22125_1','T_RS22130_1',
                             'T_RS22135_1','T_RS22140_1','STY4659_1','T_RS22150_1',
                             'T_RS22155_1', 'T_RS22160_1'),
                           c('vexE','vexD','vexC','vexB','vexA','tviE','tviD',
                             'tviC', 'tviB','tviA')))

tsl_pil <- data.frame(cbind(c('T_RS21585_1','T_RS21590_1','STY4543_1','T_RS21610_1',
                              'T_RS21615_1','T_RS21620_1','STY4547_1','SPAB_RS21825_1',
                              'T_RS21635_1'),
                            c('pilL','pilM','pilO','pilP','pilQ','pilR','pilS','pilU','pilV')))
tsl_tcf <- data.frame(cbind( c('T_RS12955_1','T_RS12950_1','T_RS12945_1','T_RS12940_1'),
                             c('tcfA','tcfB','tcfC','tcfD')))

tsl_saf <- data.frame(cbind(c('T_RS13015_1','SNSL254_RS02590_1','STMMW_03091_1','STMMW_03101_1',
                              'STMMW_03111_1','STMMW_03121_1'),
                            c('safA_typhi','safA_newport','safA_typhimurium','safB','safC','safD')))

colnames(tsl_vi) <- c("#CHROM","geneID")
colnames(tsl_pil) <- c("#CHROM","geneID")
colnames(tsl_saf) <- c("#CHROM","geneID")
colnames(tsl_tcf) <- c("#CHROM","geneID")

vi_gt_tab<- merge(tsl_vi,vi_gt_tab)
pil_gt_tab<- merge(tsl_pil,pil_gt_tab)
saf_gt_tab<- merge(tsl_saf,saf_gt_tab)
tcf_gt_tab <- merge(tsl_tcf,tcf_gt_tab)

# merge and create a table
vi_gt_tab$cluster <- "Vi_SPI7"
pil_gt_tab$cluster <- "Pil_SPI7"
saf_gt_tab$cluster <- "Saf_SPI6"
tcf_gt_tab$cluster <- "Tcf_SPI7"

full_table <- rbindlist(list(saf_gt_tab,tcf_gt_tab,vi_gt_tab,pil_gt_tab))

full_table

# identify diagnostic snps/variants?
full_table$tag <- paste0(full_table$`#CHROM`,":",full_table$POS)
full_table

# identify what loci are still variable
temp <- full_table %>% select(-c(POS,REF,`#CHROM`,cluster,geneID))
#full_table %>% filter(tag=="T_RS12945_1:737") # there was a duplicated row here- removed from initial df (ref call was CG); when the ref call is CG everyone with a C call is set to missing. 

temp <- data.frame(temp)
rownames(temp) <- temp$tag
temp$tag <- NULL
melted <- melt(t(temp))
melted_nomiss <- melted %>% filter(value != ".")
# count variants per locus
varsites <- melted_nomiss %>% group_by(Var2) %>% 
  summarise(n_unique=length(unique(value))) %>%
  filter(n_unique > 1)

varsites

full_table_varsites <- full_table %>% filter(tag %in% varsites$Var2)


m <- full_table_varsites %>% melt(id.vars = c("#CHROM",'geneID','POS','REF','cluster','tag'))
m <- m %>% filter(value != ".")
m <- m %>% mutate(val_ref_info=if_else(value==REF,"REF",value))
# add a label to each ID -> identify what clades have ref and which have alt calls
bkhd_anc <- c("IKI003","BOY008","IV3002","SJL019","SPC020","SJL020","SME005")
bkhd_mod <- c("96236","85725","PNUSAS004674","PNUSAS018445","PNUSAS007756")
typhisuis <- c("M595__NCTC20379","38K","M554__NCTC20376","SARB69CDC27768","872297")
knc139 <- c("KNC139")
ptc_mod <- c("92_0113","88_33","65_10","M422__NCTC20229","65_11bis","09_1655")
ptc_anc <- c( "TEP001","HGH1600","HGH1429","Ragna","KNL019","TAV007","GTH003","HDL004","MGL008","KNL026","KSZ005","TGP001")
choleraesuis_mod <- c("57_1_CS","RKI_16-04961","UMxr_826","UMxr_1793","0115_2009","ADRDL-14-18151","MDH-2014-00338")
choleraesuis_anc <- c("BRC034","BRC017","MGL011","ETR001","SIA006","TBS012","CPA002")
basal <- c("XBQM90","KNC095","KNC041","BPN005","KNC056","KNC181","XBQM20","CHL003")
neba <-  c("KPI001","AKY001","KKN099","GRY001","HGC040","HGC004","SUA004","NEP021","NEP018","NEP019","SAP004","KO1037","H1I004","S1L002")

unique(m$variable)
m <- m %>% mutate(lineage=if_else(variable %in% bkhd_anc, "bkhd_anc",
                                  if_else(variable %in% bkhd_mod, "bkhd_mod",
                                          if_else(variable %in% typhisuis,'typhisuis',
                                                  if_else(variable %in% knc139,'knc139',
                                                          if_else(variable %in% ptc_mod,'ptc_mod',
                                                                  if_else(variable %in% ptc_anc,'ptc_anc',
                                                                          if_else(variable %in% choleraesuis_mod,'choleraesuis_mod',
                                                                                  if_else(variable %in% choleraesuis_anc,'choleraesuis_anc',
                                                                                          if_else(variable %in% basal,'knc_xbqm_basal',
                                                                                                  if_else(variable %in% neba,'NEBA','err')))))))))))



tmp1 <- m %>% filter(val_ref_info != "REF") %>% group_by(tag) %>% 
  summarise(non_ref_ids=paste(variable,collapse=";"), 
            n_nonref=length(unique(variable)),
            nonref_lineages=paste(unique(lineage),collapse=";"),
            n_nonref_lineages=length(unique(lineage)))
tmp2 <- m %>% filter(val_ref_info == "REF") %>% group_by(tag) %>% 
  summarise(ref_ids=paste(variable,collapse=";"), 
            n_ref=length(unique(variable)),
            ref_lineages=paste(unique(lineage),collapse=";"),
            n_ref_lineages=length(unique(lineage)))
tmp3 <- merge(tmp1,tmp2)

tmp3
# save this 
#fwrite(tmp3, "variable_paraC_spi_alleles_summarised.with_lineages.csv")

# visualisation - do next NB - I think this + table could work ok here 

# x = gene order + position -> I think need to set this as a factor; y= ID; fill=genotype- make "." = "N" and set to grey.
# genotype data for visualisation 
# check - are all sites biallelic
m %>% filter(val_ref_info != 'REF') %>% select(tag,value) %>% unique() %>% group_by(tag) %>% summarise(n=n()) %>% filter(n==2)
# STMMW_03111_1:545 has two alt alleles

m %>% filter(tag=="STMMW_03111_1:545") # A is ancient bkhd only, G is all others. there is no ref calls here 
# set G to ref

m <- m %>% mutate(encoding10=if_else(tag=="STMMW_03111_1:545",if_else(value=="G",0,1),if_else(val_ref_info=="REF",0,1)))

miss <- full_table_varsites %>% melt(id.vars = c("#CHROM",'geneID','POS','REF','cluster','tag')) %>% filter(value == ".")
head(m)
head(miss)
miss <- miss %>% mutate(encoding10=".")
# need to create  data frame with chrom/gene/pos/ref/cluster/tag/variable/value/encoding01; get rid of val_ref_info and lineage
m1 <- m %>% select(-c(val_ref_info,lineage ))
m1 <- rbind(m1,miss)
# need to order tag by geneID then position, and convert to a factor with these levels

# first, need to convert geneID to a factor to be in the right order

# saf; then tcfl then pil; then vi locus
# double check which ones are in our dataset.
unique(m1$geneID) # some are "missing" - this is b/c no variation within paraC dataset, so not informative for our purposes.
# just order full set of gene IDs,
geneord <- c('safA_typhimurium','safB','safC','safD','tcfA','tcfB','tcfC','tcfD',
             'pilL','pilM','pilO','pilP','pilQ','pilR','pilS','pilU','pilV',
             'vexE','vexD','vexC','vexB','vexA','tviE','tviD','tviC', 'tviB','tviA')

m1$geneID <- factor(m1$geneID, levels=geneord)
# also need to order IDs as in tree
ordered_ids <- c("96236","85725","PNUSAS004674","PNUSAS018445","PNUSAS007756",
                 "IKI003","BOY008","IV3002","SJL019","SPC020","SJL020","SME005",
                 "S1L002","HGC040","HGC004","KPI001","AKY001","KKN099","GRY001",
                 "SUA004","NEP021","NEP018","NEP019",
                 "SAP004","KO1037","H1I004",
                 "CHL003","XBQM20", "KNC181", "KNC056","KNC041", "KNC095","BPN005","XBQM90",
                 "CPA002","TBS012","SIA006","ETR001","MGL011","BRC017","BRC034",
                 "57_1_CS","RKI_16-04961","UMxr_826","UMxr_1793","0115_2009","ADRDL-14-18151","MDH-2014-00338",
                 "KNC139","M595__NCTC20379","38K","M554__NCTC20376","SARB69CDC27768","872297",
                 "TGP001","KSZ005","KNL026","MGL008","HDL004","GTH003","TAV007","KNL019","Ragna",
                 "HGH1429","HGH1600","TEP001","92_0113","88_33","65_10","M422__NCTC20229","65_11bis","09_1655")
m1$variable <- factor(m1$variable, levels=ordered_ids)
tmp <- m1 %>% filter(variable=="TAV007")
tag_order <- tmp[order(geneID,POS,variable)]$tag
m1$tag<- factor(m1$tag,levels=tag_order)

head(m1)

ggplot(m1,aes(x=tag,y=variable,fill=encoding10)) + geom_tile() + 
  facet_wrap(~cluster, scales="free_x") + theme_cowplot()

# plot individually, this is illegible
pilin <- m1 %>% filter(cluster=="Pil_SPI7")
saf <- m1 %>% filter(cluster=="Saf_SPI6")
tcf <- m1 %>% filter(cluster=="Tcf_SPI7") # this is actually SPI-6, was a typo
vi <- m1 %>% filter(cluster=="Vi_SPI7")
ggplot(pilin,aes(x=tag,y=variable,fill=encoding10)) + geom_tile() + theme_cowplot() +
  scale_fill_manual(values=c("white","grey","black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("") + xlab("")

ggplot(vi,aes(x=tag,y=variable,fill=encoding10)) + geom_tile() + theme_cowplot() +
  scale_fill_manual(values=c("white","grey","black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("") + xlab("")

ggplot(saf,aes(x=tag,y=variable,fill=encoding10)) + geom_tile() + theme_cowplot() +
  scale_fill_manual(values=c("white","grey","black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("") + xlab("")

ggplot(tcf,aes(x=tag,y=variable,fill=encoding10)) + geom_tile() + theme_cowplot() +
  scale_fill_manual(values=c("white","grey","black")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("") + xlab("")

######### figure for pangenome alignment coverage ######### 
pg_cov_spifilt <- 'pangenome_8k.gene_coverage.for_spis.tab' # from bedtools coverage for pangenome mapping data
# filtered for SPI loci 
dat <- fread(pg_cov_spifilt)
colnames(dat)<- c("gene","start","end","n_reads_overlap","n_bases_overlap",
                  "length","pc_covered","sample_id")

id_ord_ptc <- c("EnteritidisP125109","GRY004",
                'SME005','SJL020','SJL019',"SPC020",'IV3002','IKI003','BOY008',
                'S1L002','HGC040','HGC004','KKN099','AKY001','KPI001','GRY001','H1I004',
                'SUA004','KO1037','SAP004','HO1003','NEP019','NEP021','NEP018',
                'CHL003','XBQM20','KNC181','KNC056','KNC095','KNC041','BPN005','XBQM90',
                'CPA002','SIA006','TBS012','ETR001','BRC034','BRC017',
                'MGL011','CholeraesuisSC-B67','KNC139',
                'TGP001','KSZ005','KNL026','MGL008','HDL004',
                'GTH003','TAV007','KNL019','Ragna','HGH1600',
                'Tepos_35',
                'ParatyphiC-RKS4594','TyphiTy2')
#First, get all genomic island data

gis <- fread('data/all_genomic_island_tags.tab',header=F)
colnames(gis) <- c("GI","locus_tag")
gis$gene <- paste0(gis$locus_tag,"_1")
dat_gis <- merge(dat, gis, by='gene')
length(unique(dat_gis$gene)) # 625 - some missingness (701 total)
absent_from_cap <- gis %>% filter(!gene %in% dat_gis$gene) # 75 absent from capture probes
gis_overall <- dat_gis %>% group_by(GI, sample_id) %>% summarise(total_length_considered=sum(length), total_bases_covered=sum(n_bases_overlap))
gis_overall <- gis_overall %>% mutate(pc_covered=total_bases_covered/total_length_considered)
# SPI17 not in capture probes at all
# convert to matrix for plotting
mat <- reshape2::acast(gis_overall, sample_id ~ GI, value.var = 'pc_covered')
gi_order <- c("SPI-1","SPI-2","SPI-3","SPI-4","SPI-5","SPI-6","SPI-7","SPI-8","SPI-9",
              "SPI-10","SPI-11","SPI-12","SPI-13","SPI-14","SPI-15","SPI-16","SPI-18","SPI-19",
              "SPI-20","SPI-21","CS54","SGI_1_2")
mat <- mat[rev(id_ord_ptc),gi_order]

pheatmap((mat),color=viridis::viridis(100), cluster_rows = FALSE, 
         cluster_cols = FALSE, fontsize_row = 10,fontsize_col = 10)

#SPI-3 - all lower cov
#SPI6- low cov with varn in low cov
#SPI-7 highly variable

# SPI6
spi6 <- dat_gis %>% filter(GI=="SPI-6")
# get gene level metadata including order
meta <- fread('data/spi6_locus_tags_ordered.tsv')
head(meta)
#colnames(meta)
spi6 <- merge(spi6,meta)
# get order of labels for figure
ordered_genes_spi6 <- meta[order(order)]$gene_name_figure
spi6_tmp <-  spi6 %>% select(sample_id, gene_name_figure, pc_covered) %>% unique()
spi6_mat <- reshape2::acast(spi6_tmp, sample_id ~ gene_name_figure, value.var = 'pc_covered')
spi6_mat <- spi6_mat[rev(id_ord_ptc),ordered_genes_spi6]
pheatmap((spi6_mat),color=viridis::viridis(100), cluster_rows = FALSE, 
         cluster_cols = FALSE, fontsize_row = 10,fontsize_col = 10)

# spi 7 

spi7 <- dat_gis %>% filter(GI=="SPI-7")
meta <- fread('data/spi7_locus_tags_ordered.tsv')
head(meta)

ordered_genes_spi7 <- meta[order(order_spi)]$label
spi7 <- merge(spi7,meta)
spi7_tmp <-  spi7 %>% select(sample_id, label, pc_covered) %>% unique()
spi7_mat <- reshape2::acast(spi7_tmp, sample_id ~ label, value.var = 'pc_covered')
spi7_mat <- spi7_mat[rev(id_ord_ptc),ordered_genes_spi7]
pheatmap((spi7_mat),color=viridis::viridis(100), cluster_rows = FALSE, 
         cluster_cols = FALSE, fontsize_row = 10,fontsize_col = 10)
spi7_genes <- spi7 %>% select(c(gene,colnames(meta) )) %>% unique()
spi6_genes<- spi6 %>% select(c(gene,colnames(meta) )) %>% unique()

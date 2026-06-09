rm(list=ls())
library(tidyverse)
library(data.table)
library(cowplot)


# bootstrapping analysis 
# here, subsetting to 2% of all sites 
bootmiss <- fread('bootstrap_missingness_500.tab') # output of plink --missing for 500 replicates of random subsets of sites
head(bootmiss)
# how replicable is this across boots
ggplot(bootmiss, aes(x=ID, y=f_miss)) + geom_violin() + theme_bw()
# filter for samples actually included in analysis

meta <- fread('/path/to/meta') # get sample level metadata from supplementary tables
bootmiss_qualfilt <- merge(bootmiss, meta)
ggplot(bootmiss_qualfilt, aes(x=date_mean,y=f_miss)) +
  #geom_smooth(method="lm") + 
  geom_point(aes(fill=Enterica_Tree_Position), shape=21,size=3) +  
  scale_fill_manual(values=anc_lineage_pal) + ylim(0,0.6) + theme_cowplot() + 
  xlab("Date BCE/CE") + ylab("Site missingness (bootstrapped)")

# run regression per bootstrap, not overall
models <- bootmiss_qualfilt %>% group_by(bootnum) %>% 
  do(model=lm(1-f_miss ~ date_mean, data=.))

print(models)
# Summarize coefficients by group
coefficients_summary <- models %>%
  summarise(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    pval = anova(model, test="Chisq")[1,5]
  )

ggplot(coefficients_summary, aes(y=slope,x=1)) + geom_violin()
ggplot(coefficients_summary, aes(y=pval,x=1)) + geom_violin() # all of these are like 0.1 or higher - no signif association

# model for pseudo vs mean date
anc_qualfilt <- fread('/path/to/pseudo/rates') # pseudogenisation rate table
anc_qualfilt <- merge(anc_qualfilt, meta) # to show date v pseudo rates
pseudomod <- lm(anc_qualfilter$n_pseudo/anc_qualfilter$n_genes ~ anc_qualfilter$date_mean)
summary(pseudomod)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              1.559e-02  7.196e-04  21.665  < 2e-16 ***
#   anc_qualfilter$date_mean 2.849e-06  3.573e-07   7.975 1.49e-09 ***

# how to visualise this?
ggplot(bootmiss_qualfilt, aes(x=date_mean,y=1-f_miss)) +
  #geom_smooth(method="lm") + 
  geom_point(aes(fill=Enterica_Tree_Position), shape=21,size=3) +  
  scale_fill_manual(values=anc_lineage_pal) + theme_cowplot() + 
  xlab("Date BCE/CE") + ylab("Site coverage (500 bootstraps, 2% total sites)") + ylim(0.4,1)

ggplot(coefficients_summary, aes(x=-log10(pval))) + geom_density() + theme_cowplot() + 
  geom_vline(xintercept = -log10(1.49e-09), linetype="dotted", colour="red",size=1)+ 
  theme( axis.text.y=element_blank(), axis.ticks.y=element_blank())


ggplot(coefficients_summary, aes(x=slope)) + geom_density() + theme_cowplot() + 
  geom_vline(xintercept = 2.849e-06, linetype="dotted") 

ggplot() + geom_abline(data=coefficients_summary,aes( slope=slope,intercept=intercept),alpha=0.5) + ylim(0,1) + xlim(-4000,2000)+
  geom_abline(slope=2.849e-06, intercept=1.559e-02, linetype="dotted",colour="red")

ggplot() + geom_abline(data=coefficients_summary,aes( slope=slope,intercept=1.559e-02),alpha=0.5) + ylim(0,0.1) + xlim(-4000,2000)+
  geom_abline(slope=2.849e-06, intercept=1.559e-02, linetype="dotted",colour="red",size=1) + theme_cowplot() # if you force the same intercepts



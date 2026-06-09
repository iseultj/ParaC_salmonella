# package imports
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(tidytree))
#suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(optparse))
parser <- OptionParser(description="take list of mutations annotated by snppar, and annotate whether it is fixed or not, as well as annotating number and % revertant descendants")
parser <- add_option(parser,c("-i","--infile"),action="store",help="input table")
parser <- add_option(parser,c("-t","--tr"),action="store",help="tree with internal nodes labelled")
parser <- add_option(parser,c("-o","--out"),action="store",help="output filename")
params <- parse_args(parser)

setDTthreads(1)
# annotated tree 
tr <- read.tree(params$tr)
# annotated mutations (long list)
df <- fread(params$infile,header=F)
# A464_RS09665_1 2312 N161 PNUSAS012618.fa C A anc_is_alt revertant_high
colnames(df) <- c("Gene","Pos","Anc_Node","Derived_Node","Anc","Derived","anc_classification","class")


df <- df %>% mutate(derived_interal=if_else(Derived_Node %in% tip.label(tr),"derived_tip","derived_internal"))

pseud <- df %>% filter(anc_classification=="anc_is_ref")
rev <- df %>% filter(anc_classification=="anc_is_alt")
pseud_genes <- df %>% select("Gene") %>% unique()
restab <- c("gene", "pos", "node", "is_internal", "is_fixed", "n_descendant_nodes", "revcount", "revertant_nodes","total_rev_nodes")

for (gene1 in pseud_genes$Gene) {
  sites <- pseud %>% filter(Gene==gene1) %>% select(Pos) %>% unique()
  for (pos in sites$Pos) {
    pdat <- pseud %>% filter(Gene==gene1) %>% filter(Pos==pos)
    # for each snp, get - number of pseud events (i.e. derived nodes in table)
    n_ind_events <- nrow(pdat)
    # for each of these derived nodes 
    for (node in pdat$Derived_Node) {
      # is this an internal or external node - check with number of offspring nodes
      # getting node id is nodeid(tr, node)
      descendant_nodes <- offspring(tr,nodeid(tr,node))
      if (length(descendant_nodes) == 0) {
        is_internal <- "false"
        is_fixed <- "true"
        n_descendant_nodes <- "0"
        revcount <- "NA"
        rev_vec <- "NA"
	total_rev_nodes <- "NA"
      }
      else {
        # check if any offspring nodes are in rev table
        revcount <- 0
        rev_vec <- c()
	total_rev_nodes <- 0
        for (i in descendant_nodes) {
          # need to convert node id to node label
          desclab <- nodelab(tr, i)
	n_offspring <- length(offspring(tr,i))
          n <- rev %>% filter(Gene==gene1) %>% filter(Pos==pos) %>% filter(Derived_Node==desclab) %>% nrow()
          revcount <- revcount + n
	  total_rev_nodes <- if(n > 0) total_rev_nodes + n_offspring + 1 else total_rev_nodes
          rev_vec <- if (n > 0) rbind(rev_vec, desclab) else rev_vec## this is actually just all the descendant labels
        }
        
        is_internal <- "true"
        is_fixed <- if(revcount==0) "true" else "false"
        n_descendant_nodes <- length(descendant_nodes)
      }
      result_vector <- c(gene1, pos, node, is_internal, is_fixed, n_descendant_nodes, revcount, paste(rev_vec, collapse=";"),total_rev_nodes)
      restab <- rbind(restab, result_vector)
    }
  }
}


colnames(restab) <- c("gene", "pos", "node", "is_internal", "is_fixed", "n_descendant_nodes", "revcount", "revertant_nodes","total_rev_nodes")
restab <- restab[-1,]
restab <- data.table(restab)

fwrite(restab,params$out)

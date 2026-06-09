# Steps needed to reproduce analysis for "Evolutionary history and recurrent host adaptation in ancient *Salmonella enterica*"

The following steps will allow for analyses to be reproduced for the paper "Evolutionary history and recurrent host adaptation in ancient *Salmonella enterica*" by Jackson, Neumann et al. (2025).

The following steps were run on a HPC cluster using an SGE job scheduler. Contact iseult_jackson AT eva.mpg.de for SGE batch scripts if needed, otherwise, the below commands are adaptable to other system architectures.

## 1. Data Processing: Linear Reference Genomes

Download data from ENA (Project PRJEB97698) and prepare eager input table based on library-level metadata in S. Table 7. This will also be possible to prepare with AMDirT once all metadata has been uploaded to AncientMetagenomeDir (see: <https://github.com/SPAAM-community/AncientMetagenomeDir>)

Run alignment using nf-core/eager:

```         
nextflow run nf-core/eager 
    -profile ${institutional_profile} \
    --fasta /path/to/ParatyphiC_RKS4594_NC012125.fasta \
    --bwa_index /path/to/reference_genomes/ \
    --fasta_index /path/to/ParatyphiC_RKS4594_NC012125.fasta.fai \
    --seq_dict /path/to/ParatyphiC_RKS4594_NC012125.dict  \
    --input input.tsv \
    --outdir './results_$DATE' \
    -w './work/' \
    --email $EMAIL \
    --mergedonly \
    --bwaalnl 16 \
    --bwaalnn 0.01 \
    --run_bam_filtering \
    --bam_mapping_quality_threshold 37 \
    --bam_filter_minreadlength 34 \
    --run_trim_bam \
    --bamutils_softclip \
    --bamutils_clip_single_stranded_none_udg_left 5 \
    --bamutils_clip_single_stranded_none_udg_right 5 \
    --bamutils_clip_double_stranded_none_udg_left 5 \
    --bamutils_clip_double_stranded_none_udg_right 5   \
    --bamutils_clip_single_stranded_half_udg_left 2 \
    --bamutils_clip_single_stranded_half_udg_right 2 \
    --bamutils_clip_double_stranded_half_udg_left 2 \
    --bamutils_clip_double_stranded_half_udg_right 2
```

The Paratyphi C reference genome can be replaced with the Typhi reference genome or the pSPCV plasmid genome to reproduce later analyses.

One all data has been mapped and filtered, you may want to assess authenticity, via assessing damage patterns, length distributions, edit distance distributions and heterozygosity.

Damage patterns and length distributions are assessed at a per library level: check the DamageProfiler directory for this.

For edit distance distributions, run the following:

```         
for file in $(ls *bam); do id=$(echo $file | cut -d. -f1); samtools view $file | perl -ne 'if ($_ =~ m/NM:i:(\d+)/) { print $1, "\n"}'  | sort | uniq -c | sort -k2n > ${id}.edit_dist.csv; done
```

And assess the distribution of edit distances in the output CSV using [edit_dist.R](https://github.com/LCLabTCD/aDNA_data_processing/blob/main/edit_dist.R){.uri}

Once you are satisfied with each library:

Merge BAM files to sample level using the script `merge_bamlist.py`

`python merge_bamlist.py -i input_bamlist.txt -o merged_bams`

Run duplicate removal using Picard: command format is as follows

```         
picard MarkDuplicates \
        I=${INPUTFILENAME} \
        O=${OUTPUTFILE}.duprm.bam \
        M=${OUTPUTFILE}.markdup.metrics.txt REMOVE_DUPLICATES=true
#index output of this 
samtools index  ${OUTPUTFILE}.duprm.bam
```

Finally, BAM files are softclipped (first and last 2bp of each read), and mapping quality and read length filters are applied using [filter_bam.py](https://github.com/LCLabTCD/aDNA_data_processing/blob/main/filter_bam.py)

```         
python3 filter_bam.py -i bams_to_filter.txt -mq 37 -minbp 34 -clip 2 
```

Breadth and depth of coverage of each BAM file is then checked using qualimap

`qualimap bamqc -nt 1 -bam ${INPUTFILENAME}`

Modern genome data can be split to FASTQ files of dummy reads using the script [chop_fasta_file_gzip.py](https://github.com/iseultj/killuragh_analysis_paper/blob/main/scripts/forsythia_analysis/chop_fasta_file_gzip.py) and aligned to the reference genome using nf-core/eager as above without clipping.

## 2. SNP calling

SNVs are called for phylogenetic analysis using GATK's UnifiedGenotyper and multiVCFAnalyzer.

Note that each VCF needs to go into its own directory for multiVCFAnalyzer.

A ploidy of 2 is used to allow heteroplasmic sites (for example, due to post-mortem deamination or cross-mapping), and a minimum base quality of 30 excludes softclipped bases from SNP calling.

```         
gatk2  -T UnifiedGenotyper -R /path/to/ParatyphiC_RKS4594_NC012125.fasta \
-nt 1 -glm SNP -gt_mode DISCOVERY -out_mode EMIT_ALL_SITES \
-mbq 30 -ploidy 2 -I ${INPUTFILENAME} \
--out ${OUTDIR}/${OUTPUTFILE}.ugt.mbq30.vcf.gz \
-log ${OUTDIR}/${OUTPUTFILE}.log
```

## 3. Quality Control - Heterozygosity

Each ancient VCF from step 2 should be assessed for evidence for cross-mapping and contamination by assessing heteroplasmy.

First, VCFs are parsed using `contam_ugt_vcfs.py` - this returns position, read depth, minor allele support, quality filter and mutation type (transition/transversion)

`python ./scripts/qc_scripts/contam_ugt_vcfs.py -i ${IN} -o ${OUT}.csv`

Then, summary statistics for overall heterozygosity and plots for heterozygosity in sliding windows can be calculated using `qc_summary_stats.R`. Default window size is 10kb, and step size is 1kb, but these can be changed with `-s` and `-n` respectively.

Usage:

`Rscript ./scripts/qc_scripts/qc_summary_stats.R -i ${OUT}.csv -o ${OUT} -m 3`

## 4. Multiple Sequence Alignment

Create multiple sequence alignment for modern and ancient samples passing QC filters with [multiVCFAnalyzer](https://github.com/alexherbig/MultiVCFAnalyzer) following the below command structure:

```         
java -Xmx60G -jar /path/to/MultiVCFanalyzer_0-87.jar \
        NA \ # SnpEff results file - not needed 
        /path/to/ParatyphiC_RKS4594_NC012125.fasta  \ # reference genome
        NA \ #reference genome annotation - not needed
        ./scaffold_enterica_mq37_mincov3_excl_gff_mvcfa_$DATE/ \ #output directory (needs to be created already)
        T \ #write allele frequencies - can be used to assess cross-strain mapping
        30 \ #minimum genotype quality (GATK)
        3  \ #minimum read depth supporting a call
        0.9  \ #minimal allele frequency for homozygous call
        0.9 \ #minimal allele frequency for heterozygous call 
        ParatyphiC_Regions2exclude_combined.gff \ #regions to exclude from analysis (e.g. due to high levels of cross-mapping)
        ${VCF_LIST}
```

The VCF list includes all modern and ancient VCFs to be used in phylogenetic analyses. The GFF file for masking can be found in `./data/`

## 5. Phylogenetic Analysis

Initial assessment of sample placement is performed using the multiple sequence alignment from step 4 and IQTREE for tree building.

First, the SNP alignment is filtered for maximum 5% per-site missingness using the script `./scripts/processing_scripts/Complete_partialdeletion.py`

`python3 ./scripts/processing_scripts/Complete_partialdeletion.py -f snpAlignment.fasta -o scaffold_enterica_mq37_min50pccov3x.geno05.fasta -p 95`

Then run IQTree with S. Arizonae as an outgroup to assess placement in wider S. enterica diversity:

`iqtree -s scaffold_enterica_mq37_min50pccov3x.geno05.fasta -pre scaffold_enterica_mq37_min50pccov3x.geno05.iqtree -nt 10 -mem 28G -m MFP -bb 1000 -o outgroup_Arizonae`

The above command uses 1000 rapid bootstraps to assess support, and IQTREE's ModelFinder to choose the substitution model for analysis.

Low coverage/quality ancient genomes are included in the alignment and placed in the tree one-by-one.

Most ancient genomes fall into the Para C or Birkenhead lineages (Paratyphi C, Choleraesuis, Typhisuis, Birkenhead). To investigate this more closely, multiVCFAnalyzer is run for modern and ancient genomes falling into this lineage, without the S. Arizonae outgroup, the alignment is filtered for maximum 5% per-site missingness and iqtree is run following the command structure above.

## 6. Recombination analysis

The extent to which recombination impacts the Para C/Birkenhead tree is investigated using the tool gubbins.

This tool requires the full sequence alignment (obtained from multiVCFAnalyzer). By default, this includes the reference genome used for SNP calling, which is not included in the tree, so this needs to be removed from the multiple sequence alignment for the full genome sequences before starting analysis.

```         
run_gubbins.py --starting-tree multivcfa.geno05.iqtree.treefile \
--prefix multivcfa.geno05.gubbins -c 8 \
--filter-percentage 60.0 --verbose --tree-builder iqtree \
--bootstrap 1000 --model GTRGAMMA –use-time-stamp fullAlignment.noRef.fasta
```

The tree topology and extent of recombination within this lineage can be investigated using the `plot_gubbins.R` script.

## 7. Phylotemporal Analysis

### Preparation

Steps 4-6 are re-run for modern genomes with a reported isolation date and ancient high-quality genomes (minimum 50% covered at 4X, low heteroplasmy) with radiocarbon or strong contextual dates (see Table S4 for modern genomes and S1,S2 for ancient genomes).

Full genome alignments were then masked based on the gubbins output:

```         
mask_gubbins_aln.py --aln fullAlignment.noRef.fasta \
        --gff mindp4_50pc_4xcov_c14.geno05.gubbins.recombination_predictions.gff \
        --out mindp4_50pc_4xcov_c14.masked_alignment.fasta
```

Masked and missing positions should be encoded as `N` (rather than a mix of `N` and `-`):

```         
python ./scripts/processing_scripts/fix_gaps.gubbins.py \
mindp4_50pc_4xcov_c14.masked_alignment.fasta > \
mindp4_50pc_4xcov_c14.masked_alignment.gapfix.fasta
```

SNP sites are called and monomorphic sites were counted using snp-sites:

```         
#snp call
snp-sites -m -o mindp4_50pc_4xcov_c14.masked_alignment.snps.fasta mindp4_50pc_4xcov_c14.masked_alignment.gapfix.fasta

# count
snp-sites -C mindp4_50pc_4xcov_c14.masked_alignment.gapfix.fasta > mindp4_50pc_4xcov_c14.masked_alignment.gapfix.constant_sites.txt
```

Filter SNP Alignments for maximum per-site missingness of 5%:

```         
python3 ./scripts/processing_scripts/Complete_partialdeletion.py \
-f mindp4_50pc_4xcov_c14.masked_alignment.snps.fasta \
-o mindp4_50pc_4xcov_c14.masked_alignment.snps.geno05.fasta \
-p 95
```

### Temporal signal

Root-to-tip regression is used to assess temporal signal (implemented using Treetime)

```         
treetime clock --dates ./data/metadata_for_treetime.csv --tree $iqtree_out --keep-root \
--name-column name --date_column date --plot-rtt "root_to_tip_regression.pdf" \
--aln snpAlignment.geno05.fasta
```

### BEAST setup

Uniform tip date priors for dated samples, and monophyletic constraints for the Birkenhead and Para C clades were used, and both strict and relaxed lognormal clock models were tested.

XML files to reproduce dating analysis are provided in `./data/beast_xmls/`

An example command to run analyses:

```         
for file in $(ls *xml); do prefix=$(echo $file | cut -d. -f1-2); for i in {1..2}; do seed=$(echo $RANDOM); qsub -V -b y -N beast.${prefix}.run${i} -cwd -l h_vmem=50G -pe smp 4 beast2 -seed ${seed} -packagedir /mnt/archgen/users/iseult/.beast/ -prefix ${prefix}.seed${seed}.run${i} -beagle -threads 4 ${file}; done; done
```

## 8. Pathogenicity Islands

### Presence/Absence

The presence of pathogenicity islands can be assessed by aligning to the reference genome [Typhi CT18 (NC_003198.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_003198.1), which has SPI1-10 annotated. The co-ordinates of these loci are given in BED format in Table S6.

Coverage is calculated as follows:

```         
bedtools coverage -a ${bam} -b ${spi_bed} > ${id}.SPI1-10_cov.bed
```

Specific gene co-ordinates within these pathogenicity islands can be identified using `bedtools intersect` between gene annotations for [Typhi CT18 (NC_003198.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_003198.1), and the SPI bed file, and bedtools coverage is used to calcuate per-gene coverage, as above. Bed files can be aggregated for each genome for visualisation.

This was also assessed in annotated genes from the *S. enterica* pangenome from Key et al. (2020), which is available in the `data` section of this repository along with a `bed` file with co-ordinates of SPI loci; mapping to this reference was performed slightly differently, and is outlined below. An identical `bedtools coverage` command was used to assess breadth of coverage.

### Sequence Variation

Sequence variation across key SPI loci of interest was assessed from both the Typhi reference alignment and the pangenome reference alignment.

#### Typhi Reference

Mapped BAM files should be filtered as for the Paratyphi C alignments.

Call variants against the *S. Typhi* CT18 reference genome as above, with GATK's UnifiedGenotyper with the `--EMIT_ALL_SITES` option, as for the Paratyphi C reference alignment:

```         
gatk2 -T UnifiedGenotyper \
        -R /path/to/TyphiCT18_NC003198.fasta \
        -nt 1 \
        -glm BOTH \
        -gt_mode DISCOVERY \
        -out_mode EMIT_ALL_SITES \
        -mbq 30 \
        -ploidy 2 \
        -I ${INPUTFILENAME} \
        --out ./unified_genotyper/${PREFIX}/${PREFIX}.mbq30.ugt.vcf.gz \
        -log ./unified_genotyper/${PREFIX}/${PREFIX}.mbq30.ugt.log
```

These VCFs should then be filtered for sites with a read depth within two standard deviations of the mean genomic coverage, with an absolute minimum depth of 3X and with at least 90% of reads supporting a call:

```         
python ./scripts/filter_ugt_vcfs.py -i ${input_vcf} -m ${mean_cov} -s ${std_cov} -c 3 -u 0.9 | bgzip -c > ${output_vcf}
```

Then, VCFs for individual genomes can be merged using `bcftools merge`

```         
bcftools merge --file-list merge_filtered_vcfs.txt -Oz -o modern_ancient_subsets.typhi_aln_sites.filtered.vcf.gz 
```

And filtered to regions of interest using `bcftools view`

For the pilin locus:

```         
bcftools view -S SPI7_samples.txt -r NC_003198.1:4422753-4437054 -Oz -o SPI7_pos_samples_only.pilin_locus.vcf.gz modern_ancient_subsets.typhi_aln_sites.filtered.vcf.gz
```

For the Vi capsule locus:

```         
bcftools view -S SPI7_samples.txt -r  NC_003198.1:4510582-4524679 -Oz -o SPI7_pos_samples_only.vi_locus.vcf.gz modern_ancient_subsets.typhi_aln_sites.filtered.vcf.gz
```

##### Phylogenetic Analysis

The filtered VCFs can then be used to construct a multi-FASTA file:

The script vcf2fasta can be found here: <https://github.com/iseultj/killuragh_analysis_paper/blob/main/scripts/mutans_analysis/vcf2fasta.py>

```         
# 0 - remove samples with >25% missingness and output an uncompressed VCF
# we use plink1.9 here

plink --vcf ${vcfgz_file} --mind 0.25 --recode vcf out ${vcf} 

# 1 - vcf-to-tab (from vcftools)
vcf-to-tab < ${vcf} > merged_tabbed.vcf
# vcf2fasta 
python vcf2fasta.py 25
# concatenate individual fastas -> multifasta
cat *.fa > ${prefix}.fasta
```

This multifasta is not necessarily properly aligned, particularly if there are small indels present in a subset of genomes.

Perform multiple sequence alignment on this multiFASTA file with your aligner of choice; for this manuscript we used ClustalW as implemented in SeaView (v5.0.5).

The relationship between these sequences can then be explored using the script `scripts/spis/spi_analysis.R`

##### Functional Annotations

The VCFs created by `bcftools view` (just above "Phylogenetic Analysis") can be annotated using snpEff; as these VCFs are called against the Typhi CT18 reference genome, the pre-built database for this reference genome can be used for annotations. Note that snpEff can only use an uncompressed VCF, so filtering for a minimum allele count of 1 (using bcftools) is advised.

```         
java -Xmx10G -jar /path/to/snpEff.jar eff -v -ud 100 \
-c /path/to/snpEff.config -i vcf -o vcf \
Salmonella_enterica_subsp_enterica_serovar_typhi_str_ct18 \
${input_vcf} \
> ${output_prefix}.snpEff.vcf
```

This VCF can then be inspected to assess non-synonymous variation.

### Ancestral State Reconstruction

The presence of SPI-7 at ancestral nodes in the Para C tree was assessed using ancestral state sreconstruction with `ape` and `phytools`; the script to perform this analysis can be found in `./scripts/spis/anc_state_reconstruction.seq_variation.spis.R`.

### Presence of SPIs from Pangenome Mapping

Data from the pangenome alignment (see below) can be used to assess both the presence of variable genes (which were mapped to badly in the Typhi CT18 reference), and also sequence variation across key genes in SPI-6 and -7. Mapping and variant calling steps are detailed in the "Pseudogenisation" section below. The script to create tables and figures from the resulting VCFs can be found in `./scripts/spis/anc_state_reconstruction.seq_variation.spis.R`

## 9. Pseudogenisation

### Alignment to pangenome reference

Trimmed, merged FASTQ files (produced by AdapterRemoval during sequence data processing by nf-core/eager) are used as input for mapping to the pangenome reference using `bwa mem`. For ease of processing, map libraries created using different UDG treatment separately . The flag `-L 0` is used to improve mapping to the ends of genes in the pangenome. Sorting the BAM file is performed separately due to a quirk of the system this work was performed on, but if it is possible to pipe your `samtools view` command into `samtools sort` without getting a truncated output, I would recommend doing this.

Example command:

```         
# $SM = sample ID $LB - library 
bwa mem -L 0 -R "@RG\tID:${SM}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA\tPU:dummy\tCN:dummy\tDS:ancient_udghalf\tPM:dummy" \
        -t 5 \
        ./geneID8726_cov100p_wg21k.fa \ # pangenome reference
        ${INPUTFILENAME} 2> bwa/${SM}/${OUTPUTFILE}.bwa.log \
        | samtools view -Sb -F4 -  -o  bwa/${SM}/${OUTPUTFILE}.unsorted.bam 2> bwa/${SM}/${OUTPUTFILE}.view.log
        
# sort BAM file
samtools sort -O bam -o bwa/${SM}/${OUTPUTFILE}.sorted.bam   bwa/${SM}/${OUTPUTFILE}.unsorted.bam 

# filtering: deduplication, 
picard MarkDuplicates I=${INPUTFILENAME} \
        O=./filter/${SAMPLE}/${OUTPUTFILE}.duprm.bam \
        M=./filter/${SAMPLE}/${OUTPUTFILE}.markdup.metrics.txt \
        REMOVE_DUPLICATES=true
# indel realignment, 
mkdir -p  ./indel/
#RealignerTargetCreator
gatk2 -T RealignerTargetCreator \
        -R ./geneID8726_cov100p_wg21k.fa \
        -I ${INPUTFILENAME} \
        -o ./indel/${PREFIX}.intervals
#Indel Realign
gatk2 -T IndelRealigner \
        -R ./geneID8726_cov100p_wg21k.fa \
        -I ${INPUTFILENAME} \
        -targetIntervals ./indel/${PREFIX}.intervals \
        -o ./indel/${PREFIX}.indel.bam
# mapping quality, readlength, softclipping
# n varies based on UDG treatment
# full UDG: 2; half UDG: 5; non-UDG 7

python3 filter_bam.py -i bams_to_filter.txt -mq 37 -minbp 34 -clip ${n}  

# merge bam files to sample level 

picard MergeSamFiles I=${file1} I=${file2} O=${outfile}
```

### Coverage across pangenome reference

Coverage across each locus in the pangenome reference can be assessed using `bedtools coverage`:

`bedtools coverage -a ${bam} -b data/pangenome/genome.bed > ${id}.pg_aln_cov.bed`

### Variant calling

Variants (both SNPs and indels) are called from sample-level BAM files, using UnifiedGenotyper

```         
gatk2 -T UnifiedGenotyper \
        -R ./geneID8726_cov100p_wg21k.fa \
        -nt 1 \
        -glm BOTH \
        -gt_mode DISCOVERY \
        -out_mode EMIT_ALL_SITES \
        -mbq 30 \
        -ploidy 2 \
        -I ${INPUTFILENAME} \
        --out ./02_genotyping/${OUTPUTFILE}.mbq30.ugt.vcf.gz \
        -log ./02_genotyping/${OUTPUTFILE}.mbq30.ugt.log
```

### Variant filtering

Resulting VCFs are then filtered for both minimum and maximum coverage within 2 standard deviations of the mean coverage from the Paratyphi C alignment calculated using qualimap, as well as requiring a minimum of 90% reads at a site supporting a call. Any sites not passing these filters are set to missing. This is performed using the script `filter_ugt_vcfs.py`:

```         
python filter_ugt_vcfs.py -i ${INPUTFILENAME} -m ${mean} -s ${stdev} -c 3 -u 0.9 | bgzip -c > ./filtered_vcfs/${OUTPUTFILE}.within2sd_min3.min90pc.phap.vcf.gz
```

Filtered VCFs can then be merged together using `bcftools`:

```         
bcftools merge --file-list ./complete_list_of_vcfs.excl_bad_samples.txt  --threads 5 -Oz -o ./03_input_variants/full_salmonella.full_pangenome_aln.excl_bad_samples.vcf.gz
```

The above VCF consists of all sites, which is not necessarily of interest: it should then be filtered for a minimum allele count of 1 to keep variable sites only:

```         
bcftools view -c 1 ./03_input_variants/full_salmonella.full_pangenome_aln.excl_bad_samples.vcf.gz -Oz -o ./03_input_variants/full_salmonella.full_pangenome_aln.excl_bad_samples.ac1.vcf.gz
```

At this stage, the resulting VCF can be used for the analysis of sequence variation in pathogenicity island loci (see "Presence of SPIs from Pangenome Mapping") for details.

For annotation, an additional step, atomisation, is required to split multiallelic sites across multiple lines of the VCF to make it easier to filter sites based on annotation. Additionally, an uncompressed VCF is required for snpEff, so use the flag `-Ov` to create this.

```         
bcftools norm --atomize --atom-overlaps '*' ./03_input_variants/full_salmonella.full_pangenome_aln.excl_bad_samples.ac1.vcf.gz -Ov -o ./03_input_variants/full_salmonella.full_pangenome_aln.excl_bad_samples.ac1.atomised.vcf
```

### SNPEff database building

A custom database needs to be constructed for the tool snpEff to enable variant annotation. GFF and FASTA files for the pangenome reference are provided in `./data/`.

First, create an entry in the `snpEff.config` file for the pangenome dataset, and add a codon table parameter. Example:

```         
# Salmonella pangenome covered by capture probes
salmonella_pangenome_8k_cap.genome :     salmonella_pangenome_8k_cap
    salmonella_pangenome_8k_cap.${contig}.codonTable: Bacterial_and_Plant_Plastid ## you need to add for each contig
```

Then, create a directory for this dataset, and add the GFF file to the snpEff directory:

`./snpEff/data/salmonella_pangenome_8k_cap/genes.gff`

The reference sequence is at the end of the `genes.gff` file provided, but if this were not present, the FASTA would also need to be added- for example to `./snpEff/data/genomes/salmonella_pangenome_8k_cap.fa` .

Then, the database can be built:

`java -jar snpEff.jar build -gff3 -noCheckCds -noCheckProtein -v salmonella_pangenome_8k_cap`

Once this is created, the database is ready for annotation.

### Variant annotation

The filtered, merged VCF can then be used as input to snpEff with the custom database:

```         
java -Xmx10G -jar snpEff.jar eff -v -ud 100 \
-c ./snpEff/snpEff.config -i vcf -o vcf \
salmonella_pangenome_8k_cap \
./03_input_variants/full_salmonella.full_pangenome_aln.excl_bad_samples.ac1.atomised.vcf \
> ./04_snpEff/full_salmonella.full_pangenome_aln.excl_bad_samples.ac1.atomised.snpEff.vcf
```

The annotated VCF can then be used to assess both changing pseudogenisation rates, as well as convergent pseudogenisation.

### Pseudogenisation Rates

For each individual *S. enterica* genome, pseudogenisation rates are calculated individually based on the number of genes covered ≥90% at 1X with at least one predicted pseudogenising mutation from the snpEff annotations.

First, restrict the annotated VCF to sites with predicted "HIGH" impact by snpEff: these are the predicted pseudogenising mutations.

Then, split this VCF to sample-level and restrict to sites with homozygous ALT calls (i.e. restrict to sites with these predicted pseudogenising mutations):

`bcftools +split full_dataset.pseudogenising_muts.vcf.gz -Oz -o split_vcf_dir -i 'GT="1/1"'`

Within the directory `split_vcf_dir`, there will be one VCF per sample ID, containing only predicted pseudogenising mutations.

The per-sample BED files with data on per-gene breadth of coverage (from "Coverage across pangenome reference") can then be used to aggregate pseudogenisation rates across genes classified as "present".

First, the names of all genes annotated as pseudogenes are aggregated:

```         
bcftools view ${INPUTFILENAME} | grep -v \# | awk '{print $1}' | sort | uniq -c > gene_level_aggregation/${OUTPUTFILE}.agg_gene.txt
```

This can then be filtered to genes covered \> 90% and normalised by the total number of genes covered ≥90% at ≥1X:

```         
for file in $(ls gene_level_aggregation/${OUTPUTFILE}.agg_gene.txt); do id=$(basename $file| rev | cut -d. -f3- | rev); n_genes=$(cat bedtools_cov/${id}.*bed | awk '{if ($7 > 0.9) print $1}' | wc -l); n_pseudo=$(cat bedtools_cov/${id}.*bed | awk '{if ($7 > 0.9) print $1}' | grep -f - ${file} |  wc -l); echo ${id} ${n_genes} ${n_pseudo}; done >> aggregate_pseudogene_data.pangenome_alignment.txt
```

#### Sensitivity to low-coverage data

The sensitivity of these pseudogene rate estimations to low-coverage data can be assessed using downsampling and bootstrapping analyses.

##### Downsampling

High-coverage genomes (In this study: TAV007, CPA002, KPI001) are downsampled using `samtools view –subsample` , and the full pseudogenisation rate pipeline can be run:

```         
# table of input bam, target coverage, downsample proportion as input 
input_array=($(cat input_for_downsampling_reps.csv))
inline="${input_array[$SLURM_ARRAY_TASK_ID-1]}" # move along line
id=$(echo $inline | cut -d, -f1)
bam=$(echo $inline | cut -d, -f2)
target_cov=$(echo $inline | cut -d, -f3)
ds_prop=$(echo $inline | cut -d, -f4)
rep=$(echo $inline | cut -d, -f5)

mkdir -p rep_${rep}/${target_cov}_ds/${id}/

# downsample with samtools view

# set random seed
seed=$(echo $RANDOM)
# make a record of what seed is being used
echo "${id} downsampled to ${target_cov} using seed ${seed}."

samtools view --subsample ${ds_prop} --subsample-seed ${seed} -bo rep_${rep}/${target_cov}_ds/${id}/${id}.${target_cov}_ds.pangenome_aln_filtered.bam ${bam}
```

For the downsampling analysis, identical steps are run as above to calculate pseudogenisation rates, but the VCFs aren't merged together, so the steps to split VCFs aren't necessary.

##### Bootstrapping analysis

To see if missingness might systematically reduce pseudogenisation estimates, a bootstrapping approach can be taken.

1.  Randomly subsample the VCF containing the ancient genomes used in the pseudogenisation analysis to 2% of the total variant sites in the VCF across 500 replicates.

```         
# make seperate dirs for organisation
mkdir -p boots/${bootnumb}/

# all positions must be in format chr:pos
# 2% of called sites is 10158 positions
shuf -n 10158 all_positions.txt >  boots/${bootnumb}/sites_for_boot${bootnumb}.txt

# subset to these sites


bcftools view -R  boots/${bootnumb}/sites_for_boot${bootnumb}.txt -Ov -o boots/${bootnumb}/subset.vcf full_salmonella.full_pangenome_aln.ancs_only_all.01-04-2025file.ac1.anno.vcf.gz
```

2.  Calculate individual-level missingness across these subsampled VCFs

```         
# calculate missingness

plink --vcf  boots/${bootnumb}/subset.vcf --allow-extra-chr --missing --out boots/${bootnumb}/subset.misscalc

# clean up

rm boots/${bootnumb}/subset.vcf
```

3.  Plot missingness against sample date (and perform linear regression of missingness against sample date) and compare this to the pseudogenisation rate vs sample date analysis.

```         
# create table: bootstrap index, id, n_miss, n_geno, f_miss
for file in $(ls boots/{1..500}/*imiss); do bootnum=$(echo $file | cut -d/ -f2); cat ${file} | grep -v F_MISS | awk -v b="${bootnum}" '{print b,$2,$4,$5,$6}'; done > bootstrap_missingness_500.tab
```

The script `pseudo_rates_bootstrapping.R` can then be used to analyse the results of this analysis.

### Convergent pseudogenisation

Convergent pseudogenisation is assessed using the tool SNPPar. This tool requires a tree constructed from all the samples to be analysed, as well as tables of snps to be analysed.

For the purposes of this project, we created separate tables for each gene, and split the table by snpEff annotation (High impact (i.e. pseudogenising), missense, and synonymous).

First, a maximum likelihood tree can be constructed as above (5. Phylogenetic Analysis) from the SNP alignment generated from mapping reads to the Paratyphi C reference. Make sure that the sample list for the SNPPar analysis and the tree match exactly, that bootstrap support values are removed, and that the sample IDs are identical.

Second, the set of genes to perform the analysis on needs to be defined - take those covered ≥90% (from `bedtools coverage`) in ≥95% of the samples to be used for analysis.

Third, genbank files for each gene is required: this can be created from the `genes.gff` file, for example using the following perl script: <https://gist.github.com/avrilcoghlan/5386810> .

Fourth, the input tables for SNPPar need to be created from the input VCFs. The below code is for an SGE array of genes to analyse:

```         
while read line; do input_array+=(${line}); done < ./s_enterica_core_genes.defined_by_dataset.10-04-2025.txt 

# loop through array

scaffold="${input_array[SGE_TASK_ID -1]}"
# prep directories

mkdir -p ./05_snppar_prep/split_vcfs/high
mkdir -p ./05_snppar_prep/split_vcfs/missense
mkdir -p ./05_snppar_prep/split_vcfs/synonymous

mkdir -p ./05_snppar_prep/intermediate_tab/high
mkdir -p ./05_snppar_prep/intermediate_tab/missense
mkdir -p ./05_snppar_prep/intermediate_tab/synonymous

mkdir -p ./05_snppar_prep/snptable/high
mkdir -p ./05_snppar_prep/snptable/missense
mkdir -p ./05_snppar_prep/snptable/synonymous

## step 1: subset input VCF.GZs to scaffold level VCF

bcftools view ./05_snppar_prep/pangenome_aln.atomised.snpEff.subset.ids.AC1.HIGH.vcf.gz \
${scaffold} -Ov -o ./05_snppar_prep/split_vcfs/high/${scaffold}.vcf

bcftools view ./05_snppar_prep/pangenome_aln.atomised.snpEff.subset.ids.AC1.missense.vcf.gz \
${scaffold} -Ov -o ./05_snppar_prep/split_vcfs/missense/${scaffold}.vcf

bcftools view  ./05_snppar_prep/pangenome_aln.atomised.snpEff.subset.ids.AC1.synonymous.vcf.gz \
${scaffold} -Ov -o ./05_snppar_prep/split_vcfs/synonymous/${scaffold}.vcf

## step 2: run vcf-to-tab - part of vcftools
vcf-to-tab < ./05_snppar_prep/split_vcfs/high/${scaffold}.vcf > ./05_snppar_prep/intermediate_tab/high/${scaffold}.tab

vcf-to-tab < ./05_snppar_prep/split_vcfs/missense/${scaffold}.vcf > ./05_snppar_prep/intermediate_tab/missense/${scaffold}.tab

vcf-to-tab < ./05_snppar_prep/split_vcfs/synonymous/${scaffold}.vcf > ./05_snppar_prep/intermediate_tab/synonymous/${scaffold}.tab

## step 3: manipulate tabbed vcf to input snp table for snppar
## replace ^I with tab. everything else incl regexp should stay the same. 
cat ./05_snppar_prep/intermediate_tab/high/${scaffold}.tab | sed 's/\.\//N\//g' | sed 's/\/\./\/N/g' | sed 's/\//}/g' | sed 's/}\w\+^I/^I/g' | sed 's/}^I/^I/g' | sed 's/,/_/g' | sed 's/}\w\+$//g' | sed 's/^I/,/g' | sed 's/ /,/g' | sed 's/,N,/,-,/g' | sed 's/,N,/,-,/g' | sed 's/,N$/,-/g' | cut -d, -f2,4- | sed 's/POS/Pos/g' > ./05_snppar_prep/snptable/high/${scaffold}.high.snpTable.csv

cat ./05_snppar_prep/intermediate_tab/missense/${scaffold}.tab |  sed 's/\.\//N\//g' | sed 's/\/\./\/N/g' | sed 's/\//}/g' | sed 's/}\w\+^I/^I/g' | sed 's/}^I/^I/g' | sed 's/,/_/g' | sed 's/}\w\+$//g' | sed 's/^I/,/g' | sed 's/ /,/g'  | sed 's/,N,/,-,/g' | sed 's/,N,/,-,/g' | sed 's/,N$/,-/g' |   cut -d, -f2,4- | sed 's/POS/Pos/g' > ./05_snppar_prep/snptable/missense/${scaffold}.missense.snpTable.csv

cat ./05_snppar_prep/intermediate_tab/synonymous/${scaffold}.tab |  sed 's/\.\//N\//g' | sed 's/\/\./\/N/g' | sed 's/\//}/g' | sed 's/}\w\+^I/^I/g' | sed 's/}^I/^I/g' | sed 's/,/_/g' | sed 's/}\w\+$//g' | sed 's/^I/,/g' | sed 's/ /,/g'  | sed 's/,N,/,-,/g' | sed 's/,N,/,-,/g' | sed 's/,N$/,-/g' |   cut -d, -f2,4- | sed 's/POS/Pos/g' > ./05_snppar_prep/snptable/synonymous/${scaffold}.synonymous.snpTable.csv

## step 4: clean up files that need to be compressed

bgzip ./05_snppar_prep/split_vcfs/high/${scaffold}.vcf  
bgzip ./05_snppar_prep/intermediate_tab/high/${scaffold}.tab

bgzip ./05_snppar_prep/split_vcfs/missense/${scaffold}.vcf  
bgzip ./05_snppar_prep/intermediate_tab/missense/${scaffold}.tab

bgzip ./05_snppar_prep/split_vcfs/synonymous/${scaffold}.vcf  
bgzip ./05_snppar_prep/intermediate_tab/synonymous/${scaffold}.tab
```

Then, SNPPar can be run on these tables; as with the table preparation, an array is used to store the core genes to run this analysis on.

```         
while read line; do input_array+=(${line}); done < ./s_enterica_core_genes.defined_by_dataset.10-04-2025.txt 
scaffold="${input_array[SGE_TASK_ID -1]}"
mkdir -p ./06_snppar/high/${scaffold}
mkdir -p ./06_snppar/missense/${scaffold}
mkdir -p ./06_snppar/synonymous/${scaffold}

#high
snppar -s ./05_snppar_prep/snptable/high/${scaffold}.high.snpTable.csv -E simple \
        -t rooted_enterica_tree.for_snppar.no_bootstraps.nwk \
        -g ./genbank_files/${scaffold}.gbk \
        -d ./06_snppar/high/${scaffold} \
        -p ${scaffold} \
        -A 

# missense

snppar -s ./05_snppar_prep/snptable/missense/${scaffold}.missense.snpTable.csv -E simple \
        -t rooted_enterica_tree.for_snppar.no_bootstraps.nwk  \
        -g ./genbank_files/${scaffold}.gbk \
        -d ./06_snppar/missense/${scaffold}  \
        -p ${scaffold} \
        -A 

# synonymous

snppar -s ./05_snppar_prep/snptable/synonymous/${scaffold}.synonymous.snpTable.csv -E simple \
        -t rooted_enterica_tree.for_snppar.no_bootstraps.nwk  \
        -g ./genbank_files/${scaffold}.gbk \
        -d ./06_snppar/synonymous/${scaffold} \
        -p ${scaffold} \
        -A 
```

The output of these jobs can then be aggregated to CSVs, and each analyzed mutation can be annotated as being revertant or not (i.e. whether the ancestral state is the reference or alternate call) - this is important for understanding which annotated mutation events are pseudogenising or not.

```         
# aggregation by impact
for file in $(ls high/*/*all_mutation_events.tsv); do gene=$(echo $file | cut -d/ -f2); cat ${file} | grep -v Position | awk -v gene="${gene}" '{print gene,$1,$3,$4,$5,$6}'; done >> snppar_HIGH_impact.all_variants_aggregated.txt
for file in $(ls missense/*/*all_mutation_events.tsv); do gene=$(echo $file | cut -d/ -f2); cat ${file} | grep -v Position | awk -v gene="${gene}" '{print gene,$1,$3,$4,$5,$6}'; done >> snppar_missense_impact.all_variants_aggregated.txt
for file in $(ls synonymous/*/*all_mutation_events.tsv); do gene=$(echo $file | cut -d/ -f2); cat ${file} | grep -v Position | awk -v gene="${gene}" '{print gene,$1,$3,$4,$5,$6}'; done >> snppar_synonymous_impact.all_variants_aggregated.txt 

## aggregate which alleles are ref/alt per snp (split by annotated impact)
zcat intermediate_tab/high/*tab.gz  | grep -v \# | awk '{print $1,$2,$3}' > high_impact_snps.ref_info.txt
zcat intermediate_tab/missense/*tab.gz  | grep -v \# | awk '{print $1,$2,$3}' > missense_impact_snps.ref_info.txt
zcat intermediate_tab/synonymous/*tab.gz  | grep -v \# | awk '{print $1,$2,$3}' > synonymous_impact_snps.ref_info.txt

## check whether the ancestral state is ref or alt
while read line; do gene=$(echo $line | awk '{print $1}'); pos=$(echo $line | awk '{print $2}'); anc=$(echo $line | awk '{print $5}'); ref=$(cat high_impact_snps.ref_info.txt | grep ${gene} | grep ${pos} | awk '{print $3}'); if [[ ${anc} = ${ref} ]]; then echo ${line} anc_is_ref high; else echo ${line} anc_is_alt revertant_high; fi; done < snppar_HIGH_impact.all_variants_aggregated.txt >> full_info_from_snppar.high_impact.txt ## and similarly for missense and synonymous sites. 

## annotate % revertant descendant nodes from pseudogenising mutations
# the input tree for this script is any of the output trees from SNPPar with internal node labels
Rscript annotate_fixed.R -i full_info_from_snppar.high_impact.txt  -t input_tree.labelled_nodes.nwk -o snppar_high_impact.annotated.tab
```

This annotated table can then be used as input for the R scripts used to perform convergent pseudogenisation analysis.

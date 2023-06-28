library(data.table)
library(CLONETv2)
library(TPES)
library(ggplot2)


normal = fread("Control.csv",data.table=F)
normal$af = normal$altCount/normal$totalCount
tumor = fread("Tumor.csv",data.table=F)
tumor$af = tumor$altCount/tumor$totalCount

pileup.normal = normal[,c(1,2,4,5,14,8)]
colnames(pileup.normal) = c("chr","pos","ref","alt","af","cov")

pileup.tumor = tumor[,c(1,2,4,5,14,8)]
colnames(pileup.tumor) = c("chr","pos","ref","alt","af","cov")

seg.tb <- fread("~/Desktop/HumanGenomics/Project_2023/copy_number/SCNA.copynumber.called.seg",data.table=F)

setwd("~/Desktop/HumanGenomics/Project_2023/purity_ploidy")

bt <- compute_beta_table(seg.tb, pileup.tumor, pileup.normal)

## Compute ploidy table with default parameters
pl.table <- compute_ploidy(bt)

adm.table <- compute_dna_admixture(beta_table = bt, ploidy_table = pl.table)

allele_specific_cna_table <- compute_allele_specific_scna_table(beta_table = bt,
                                                                ploidy_table = pl.table, 
                                                                admixture_table = adm.table)


check.plot <- check_ploidy_and_admixture(beta_table = bt, ploidy_table = pl.table,
                                         admixture_table = adm.table)
print(check.plot)

ggplot(data=allele_specific_cna_table,aes(x=cnA,y=cnB,col=log2.corr)) +
coord_fixed() +
geom_point()


# TPES
snv.reads = fread("../somatic_variant/somatic.pm",data.table=F)
snv.reads = snv.reads[which(snv.reads$somatic_status=="Somatic"),]
snv.reads = snv.reads[,c("chrom","position","position","tumor_reads1","tumor_reads2")]
colnames(snv.reads) = c("chr","start","end","ref.count","alt.count")
snv.reads$sample = "Sample.1"

TPES_purity(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)

TPES_report(ID = "Sample.1", SEGfile = seg.tb,
            SNVsReadCountsFile = snv.reads,
            ploidy = pl.table,
            RMB = 0.47, maxAF = 0.6, minCov = 10, minAltReads = 10, minSNVs = 1)


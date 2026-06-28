library(data.table)
library(dplyr)


#------------- LD clumping ---------------------

phen = 'acc'
gwas = fread("fastGWA.fastGWA")

library(ieugwasr)
bfile = "/data2/xsq/1kg.v3/EUR"  # reference file for LD computation

gwas = gwas[gwas$P<0.05,]
gwas = dplyr::rename(gwas, rsid=SNP, pval=P)
gwas$SNP = paste0(gwas$rsid, '_', gwas$A1,gwas$A2)

gwas_clump = ld_clump(dat = gwas, 
                    clump_kb = 5000,
                    clump_r2 = 0.1,
                    clump_p = 1,
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    bfile = bfile,
                    pop = "EUR"
)




#------------- PRS ---------------------

cutoff= 0.05

prs_sum=numeric()
for(chr in 1:22){
  # chr=5
  filename = paste0("prs/",phen,"/clump_",cutoff,"/chr",
                    chr,".sscore")
  if(!file.exists(filename)){next}
  prs <- read.delim(filename)
  prs_sum = cbind(prs_sum, prs$SCORE1_SUM)
}
ncol(prs_sum)

prs_sum = cbind(FID=prs$X.FID, IID=prs$IID, PRS=rowSums(prs_sum))





















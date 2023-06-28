library(EthSEQ)

setwd("G:/Il mio Drive/EthSEQ")

ethseq.Analysis(
  target.vcf = "./Control.VARSCAN.vcf",
  out.dir = "./Control",
  model.gds = "SS2.Light.Model.gds",
  cores=1,
  mbq = 20,
  mrq = 1,
  mdc = 10,
  run.genotype = TRUE,
  verbose=TRUE,
  composite.model.call.rate = 1,
  space="3D")

ethseq.Analysis(
  target.vcf = "./Tumor.VARSCAN.vcf",
  out.dir = "./Tumor",
  model.gds = "SS2.Light.Model.gds",
  cores=1,
  mbq = 20,
  mrq = 1,
  mdc = 10,
  run.genotype = TRUE,
  verbose=TRUE,
  composite.model.call.rate = 1,
  space="3D")


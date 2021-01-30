

  library(divergence.preSE)
  library(rutils)
  library(GSEABase)

  data.pheno.old = read.table(gzfile("DATA/TCGA_BRCA_Pheno.csv.gz"),
                              sep=",", header=TRUE)
  rownames(data.pheno.old) = data.pheno.old$sample
  
  data.pheno = read.csv("DATA/TCGA_BRCA_clinicalMatrix.gz", 
                        sep="\t", as.is=TRUE)  
  rownames(data.pheno) = gsub("-", ".", as.character(data.pheno$sampleID))
  
  load("DATA/geneSetCollectionSet.MSigDBv6.1.rda")
  genesets = lapply(geneSetCollectionSet, geneIds)
  
  # process expression data
  
  data.exp = read.csv("DATA/TCGA_BRCA_HiSeqV2.gz", sep="\t")
  
  mat.exp = data.matrix(data.exp[, -1])
  rownames(mat.exp) = as.character(data.exp[, 1])
  
  sel.samples = intersect(rownames(data.pheno), colnames(mat.exp))
  
  mat.exp = mat.exp[, sel.samples]
  pheno.exp = data.pheno[sel.samples, ]
  
  sel.normal = which(pheno.exp$sample_type == "Solid Tissue Normal")
  
  sel.tumor = which(pheno.exp$sample_type == "Primary Tumor")

  div.exp = computeUnivariateDigitization(Mat=mat.exp[, sel.tumor], 
                                      baseMat=mat.exp[, sel.normal], 
                                      computeQuantiles=TRUE, 
                                      gamma=1:5/10)
  
  genesets.KEGG = lapply(genesets$Broad.c2.CP.KEGG, function(x) intersect(x, rownames(mat.exp)))
  
  div.exp.KEGG = computeMultivariateDigitization(Mat=mat.exp[, sel.tumor], 
                                                 baseMat=mat.exp[, sel.normal], 
                                                 FeatureSets=genesets.KEGG,
                                                 computeQuantiles=TRUE, 
                                                 gamma=1:5/10)
  
  genesets.CGP = lapply(genesets$Broad.c2.CGP, function(x) intersect(x, rownames(mat.exp)))
  
  div.exp.CGP = computeMultivariateDigitization(Mat=mat.exp[, sel.tumor], 
                                                 baseMat=mat.exp[, sel.normal], 
                                                 FeatureSets=genesets.CGP,
                                                 computeQuantiles=TRUE, 
                                                 gamma=1:5/10)
  
  
  save(data.pheno, data.pheno.old, mat.exp, pheno.exp, div.exp,
       genesets, genesets.KEGG, genesets.CGP,
       div.exp.KEGG, div.exp.CGP,
       file="data_1.rda")
  
  # process methylation data
  
  data.met = read.csv("DATA/TCGA_BRCA_HumanMethylation450.gz", sep="\t")
  
  mat.met = data.matrix(data.met[, -1])
  rownames(mat.met) = as.character(data.met[, 1])

  # remove NAs
  r = which(apply(mat.met, 1, function(x) any(is.na(x))))
  mat.met = mat.met[-r, ]
  
  # all(is.finite(mat.met))
  
  sel.samples = intersect(rownames(data.pheno), colnames(mat.met))
  
  mat.met = mat.met[, sel.samples]
  pheno.met = data.pheno[sel.samples, ]
  
  sel.normal = which(pheno.met$sample_type == "Solid Tissue Normal")
  
  sel.tumor = which(pheno.met$sample_type == "Primary Tumor")
  
  div.met = computeUnivariateDigitization(Mat=mat.met[, sel.tumor], 
                                          baseMat=mat.met[, sel.normal], 
                                          computeQuantiles=TRUE, 
                                          gamma=1:9/10)
  
  # cpg-gene annotations
  
  data.cpg.ann = read.csv("DATA/CpGsAnnotatedToGenesHg19.txt.gz", 
               sep="\t", as.is=TRUE)
  
  genes.ann = unique(data.cpg.ann$annot.symbol)
  
  genes.cpg.map = lapply(genes.ann, function(x) unique(data.cpg.ann$id[which(data.cpg.ann$annot.symbol == x)]) )
  names(genes.cpg.map) = genes.ann
  
  genes.cpg.map.2 = lapply(genes.cpg.map, function(x) intersect(x, rownames(mat.met)))
  
  genes.cpg.map.2_uni = genes.cpg.map.2[sapply(genes.cpg.map.2, length) == 1]
  genes.cpg.map.2_multi = genes.cpg.map.2[sapply(genes.cpg.map.2, length) > 1]
  
  div.met.gene = computeMultivariateDigitization(Mat=mat.met[, sel.tumor], 
                                               baseMat=mat.met[, sel.normal], 
                                               FeatureSets=genes.cpg.map.2_multi,
                                               computeQuantiles=TRUE, 
                                               gamma=c(1, 5, 7, 9)/10)

  
  cpgs.all = unique(data.cpg.ann$id)
  
  cpg.ann.list = lapply(cpgs.all, function(x){
    y = which(data.cpg.ann$id == x)
    list(
      cpg = x,
      chr = data.cpg.ann$seqnames[ y[1] ],
      start = data.cpg.ann$start[ y[1] ],
      annot.strand = unique(data.cpg.ann$annot.strand[y]),
      annot.width = unique(data.cpg.ann$annot.width[y]),
      annot.symbol = unique(data.cpg.ann$annot.symbol[y]),
      annot.type = unique(data.cpg.ann$annot.type[y]),
      annot.dist = unique(data.cpg.ann$DistanceFromAnn[y])
    )
  })
  names(cpg.ann.list) = cpgs.all
  
  cpg.genes.all = unique(Reduce(c, lapply(cpg.ann.list, function(x) x$annot.symbol)))
  
  save(div.met, data.cpg.ann, genes.cpg.map, 
       genes.cpg.map.2, genes.cpg.map.2_uni, genes.cpg.map.2_multi,
       div.met.gene,
       cpg.ann.list, cpg.genes.all,
       file="data_2.rda")
  
  save(mat.met, sel.tumor, sel.normal, genes.cpg.map.2_multi, file="data_3.rda")
  









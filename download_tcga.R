
  library(TCGAbiolinks)

  # create folders
  dir.create("BRCA_EXP")
    
  dir.create("BRCA_MET")
  
  # download expression data
  
  query.exp = GDCquery(project = "TCGA-BRCA",
                        legacy = TRUE,
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        experimental.strategy = "RNA-Seq")
  
  GDCdownload(query.exp,
              directory = "BRCA_EXP")
  
  
  
  # download methylation data
  
  query.met = GDCquery(project = "TCGA-BRCA",
                        legacy = TRUE,
                        data.category = "DNA methylation",
                        platform = "Illumina Human Methylation 450")
  
  GDCdownload(query.met,
              directory = "BRCA_MET")
  
  # process and save summarized experiment
  
  data.exp = GDCprepare(query = query.exp, directory = "BRCA_EXP")
  
  save(data.exp, file="data.exp.rda")
  
  data.met = GDCprepare(query = query.met, directory = "BRCA_MET")
  
  save(data.met, file="data.met.rda")
  
  
  # ==================
  
  library(SummarizedExperiment)

  mat.exp = log2(assays(data.exp)[[1]]+1)
  pheno.data = as.data.frame(colData(data.exp))
  
  save(mat.exp, pheno.data, file="DATA/biolinks_BRCA_exp.rda")
  
  
  mat.met = assays(data.met)[[1]]
  pheno.data = as.data.frame(colData(data.met))
  
  save(mat.met, pheno.data, file="DATA/biolinks_BRCA_met.rda")

  
  
    
  






  library(divergence.preSE)
  library(GSEABase)

  M1 = data.matrix(read.csv("DATA/GPL96_TUMOR_Expression.csv.gz"))
  P1 = read.csv("DATA/GPL96_TUMOR_Pheno.csv.gz")
  
  rownames(P1) = P1$sample
  
  l1 = load("data_1.rda")

  genes = rownames(M1)
   
  genesets.KEGG = lapply(genesets$Broad.c2.CP.KEGG, function(x) intersect(x, genes) )

  sel.n = which(P1$PAM50 == "Normal")
  
  sel.t = which(P1$PAM50 %in% c("Luminal_A", "Luminal_B"))
  groups.mic.PAM50 = P1$PAM50[sel.t]
  names(groups.mic.PAM50) = rownames(P1)[sel.t]
  
  groups.mic.PAM50 = gsub("_", " ", groups.mic.PAM50)
  
  groups.all = P1$PAM50[-sel.n]
  names(groups.all) = rownames(P1)[-sel.n]
  
  groups.all = gsub("_", " ", groups.all)
  
  gammas = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9)
  
  M1Q = computeQuantileMatrix(M1)
  
  sel.other = which(! groups.all %in% c("Luminal A", "Luminal B"))
  
  divp.exp = rowMeans(div.exp.KEGG$Mat.div)
  
  div.cor = data.frame(t(sapply(gammas, function(g){
    
    div.temp = computeMultivariateDigitization(Mat=M1Q[genes, names(groups.all)[sel.other]], 
                                        baseMat=M1Q[genes, sel.n], 
                                        FeatureSets=genesets.KEGG,
                                        computeQuantiles=FALSE, 
                                        beta=0.75,
                                        gamma=g)
    
    divp.mic = rowMeans(div.temp$Mat.div)
    
    r = cor(divp.exp[intersect(names(divp.exp), names(divp.mic))], divp.mic[intersect(names(divp.exp), names(divp.mic))])
    
    c(gamma=g, r=r)
    
  })))
  
  sel.g = div.cor$gamma[which.max(div.cor$r)]
  
  div.mic.KEGG = computeMultivariateDigitization(Mat=M1Q[genes, names(groups.all)], 
                                                 baseMat=M1Q[genes, sel.n], 
                                                 FeatureSets=genesets.KEGG,
                                                 computeQuantiles=FALSE, 
                                                 beta=0.75,
                                                 gamma=sel.g)
  
  
  
  save(div.mic.KEGG, groups.mic.PAM50, groups.all, file="data_4.rda")
  
  









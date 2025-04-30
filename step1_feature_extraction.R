# install.packages("protr")
library(protr)

rm(list=ls())


  proSeq <- readFASTA("TF_sequences.fasta")
  proSeq<- proSeq[(sapply(proSeq, protcheck))]##进行氨基酸类型完整性检查并删除非标准序列
  number=length(proSeq)

  res1 <- t(sapply(proSeq, extractProtFP, index=c(1:531),pc = 10,lag=5, silent = FALSE))
  ##res包含lag*p^2个特征数据,包括pc个主成分的自相关和每个主成分与其它主成分之间的交叉相关。已验证过。
  res2 <- t(sapply(proSeq, extractDescScales,propmat = "AATopo", index=NULL,pc = 10,lag=5, silent = FALSE))
  ##index=NULL表示选择所有特征参数。
  res3 <- t(sapply(proSeq, extractDescScales,propmat = "AAMolProp", index=NULL,pc = 10,lag=5, silent = FALSE))
  res4 <- t(sapply(proSeq, extractDescScales,propmat = "AAGeom", index=NULL,pc = 10,lag=5, silent = FALSE))
  res5 <- t(sapply(proSeq, extractDescScales,propmat = "AAConn", index=NULL,pc = 10,lag=5, silent = FALSE))
  res6 <- t(sapply(proSeq, extractDescScales,propmat = "AAEdgeAdj", index=NULL,pc = 10,lag=5, silent = FALSE))
  res7 <- t(sapply(proSeq, extractDescScales,propmat = "AAInfo", index=NULL,pc = 10,lag=5, silent = FALSE))
  res8 <- t(sapply(proSeq, extractDescScales,propmat = "AAMOE2D", index=NULL,pc = 10,lag=5, silent = FALSE))
  res9 <- t(sapply(proSeq, extractDescScales,propmat = "AAMOE3D", index=NULL,pc = 10,lag=5, silent = FALSE))
  res10 <- t(sapply(proSeq, extractDescScales,propmat = "AAWalk", index=NULL,pc = 10,lag=5, silent = FALSE))
  res11 <- t(sapply(proSeq, extractDescScales,propmat = "AAWHIM", index=NULL,pc = 10,lag=5, silent = FALSE))
  res12 <- t(sapply(proSeq, extractDescScales,propmat = "AAACF", index=NULL,pc = 10,lag=5, silent = FALSE))
  res13 <- t(sapply(proSeq, extractDescScales,propmat = "AABurden", index=NULL,pc = 10,lag=5, silent = FALSE))
  res14 <- t(sapply(proSeq, extractDescScales,propmat = "AACPSA", index=NULL,pc = 10,lag=5, silent = FALSE))
  res15 <- t(sapply(proSeq, extractDescScales,propmat = "AADescAll", index=NULL,pc = 10,lag=5, silent = FALSE))
  res16 <- t(sapply(proSeq, extractDescScales,propmat = "AA2DACOR", index=NULL,pc = 10,lag=5, silent = FALSE))
  res17 <- t(sapply(proSeq, extractDescScales,propmat = "AA3DMoRSE", index=NULL,pc = 10,lag=5, silent = FALSE))
  res18 <- t(sapply(proSeq, extractDescScales,propmat = "AAConst", index=NULL,pc = 10,lag=5, silent = FALSE))
  res19 <- t(sapply(proSeq, extractDescScales,propmat = "AAEigIdx", index=NULL,pc = 10,lag=5, silent = FALSE))
  res20 <- t(sapply(proSeq, extractDescScales,propmat = "AAFGC", index=NULL,pc = 10,lag=5, silent = FALSE))
  res21 <- t(sapply(proSeq, extractDescScales,propmat = "AAGETAWAY", index=NULL,pc = 10,lag=5, silent = FALSE))
  res22 <- t(sapply(proSeq, extractDescScales,propmat = "AARandic", index=NULL,pc = 10,lag=5, silent = FALSE))
  res23 <- t(sapply(proSeq, extractDescScales,propmat = "AARDF", index=NULL,pc = 10,lag=5, silent = FALSE))
  res24 <- t(sapply(proSeq, extractDescScales,propmat = "AATopoChg", index=NULL,pc = 10,lag=5, silent = FALSE))

  
  x1 <- t(sapply(proSeq,extractAAC))
  x2 <- t(sapply(proSeq,extractDC))
  x3 <- t(sapply(proSeq,extractTC))
  x4 <- t(sapply(proSeq,extractMoreauBroto))
  x5 <- t(sapply(proSeq,extractMoran))
  x6 <- t(sapply(proSeq,extractGeary))
  x7 <- t(sapply(proSeq,extractCTDC))
  x8 <- t(sapply(proSeq,extractCTDT))
  x9 <- t(sapply(proSeq,extractCTDD))
  x10 <- t(sapply(proSeq,extractCTriad))
  x11 <- t(sapply(proSeq,extractSOCN))
  x12 <- t(sapply(proSeq,extractQSO))
  x13 <- t(sapply(proSeq,extractPAAC))
  x14 <- t(sapply(proSeq,extractAPAAC))
  df.integrated <- cbind(res1,res2,res3,res4,res5,res6,res7,res8,res9,res10,res11,res12,res13,res14,res15,res16,res17,res18,res19,res20,res21,res22,res23,res24,x1,x2,x3,x4,x5,x6,x7,x8,x9,x11,x12,x13,x14)
  write.table(df.integrated,"all_feature_protein_completeSeq.txt",quote=F,sep ="\t",row.names = T,col.names=T)





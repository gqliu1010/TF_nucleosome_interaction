
rm(list=ls())
combined_data <- read.table("all_feature_protein_completeSeq.txt",sep ="\t",row.names = 1,header = T)

filePath <- list.files(pattern = "integrated_coef*", full.names = TRUE)
filePath <- filePath[c(2,3,5,7)]
# 使用lapply读取所有文件的内容，并合并到一个列表中
data_list <- lapply(filePath, read.table, header = F)
# 使用do.call和cbind按列追加合并
combined_featureName <- unlist(do.call(rbind, data_list))
##do.call 允许你将一个函数应用到一个列表的元素上。它的基本用法是将一个函数和一个列表作为输入，然后将列表中的元素作为参数传递给函数。这在处理多个对象时非常方便，尤其是当需要将多个对象合并或应用某个函数时。
combined_featureName <- unique(combined_featureName)
selected_feature <- combined_data[,combined_featureName]

write.table(combined_featureName, "selected_featureName_protein_completeSeq.txt",quote=F,sep ="\t",row.names = T,col.names = T)
write.table(selected_feature, "selected_feature_protein_completeSeq.txt",quote=F,sep ="\t",row.names = T,col.names = T)

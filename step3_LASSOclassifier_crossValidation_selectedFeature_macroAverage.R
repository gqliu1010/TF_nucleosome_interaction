
  library(caret) 
  library(pROC) 
  library(e1071)
  library(openxlsx)
  library(glmnet)
  library(PRROC)

  
  
  rm(list=ls())

  bindMode <- c("End.preference","Periodic.preference","Groove.preference","Dyad.preference","Gyre.spanning","Orientational.preference","Nucleosome.stability")
  raw=read.table("selected_feature_protein_completeSeq.txt",header =T,sep='\t',row.names=1,stringsAsFactors=F, comment.char = "!")
  raw$nucleosome_binding_status=NULL

  group_file=read.table("multi.label.txt",header =T,sep='\t',row.names=1,stringsAsFactors=F, comment.char = "!")
  group_file=group_file[match(rownames(raw),rownames(group_file)),]
  
  group_file$Groove.preference[group_file$Groove.preference=="Major"]="Yes"
  group_file$Groove.preference[group_file$Groove.preference=="Minor"]="No"
  group_file$Nucleosome.stability[group_file$Nucleosome.stability=="Stabilizer"]="Yes"
  group_file$Nucleosome.stability[group_file$Nucleosome.stability=="Destabilizer"]="No"
  
  
  
  cols_with_na <- colSums(is.na(raw)) > 0
  #names(data)[cols_with_na]
  raw <- raw[,!cols_with_na]
  
  for(k in bindMode){
    print(k)
    group_inf=group_file[which(group_file[,k]=="Yes" | group_file[,k]=="No"),]
    data=raw[which(group_file[,k]=="Yes" | group_file[,k]=="No"),]
    if((all(rownames(data)==rownames(group_inf)))==FALSE){
      stop("数据顺序有错误")
    }
    
    data$tag <- ifelse(group_inf[,k]=="Yes",1,0)
    #if(length(data$tag[data$tag==1])<20){next}
    set.seed(2)
    M=10 #M-fold cross validation
    # fraction=1/M
    # ind<-sample(M,nrow(data),replace = T,prob = c(rep(fraction,M)))
    ########正负样本混起来一起分组，当正样本（或负样本）较少时，这种分组方法可能会导致划分的测试集或测试集中没有正样本（或负样本）############################################################
    ##folds <- createFolds(data$tag, k = 10) # 创建10折交叉验证的分层抽样
    
    #######正负样本分别分组，以保证数量较少的正或负样本在每次的测试集和训练集中至少出现1例。
    posIND=which(data$tag==1)
    negIND=which(data$tag==0)
    ##library(caret)
    posfolds <- createFolds(posIND, k = M,list=T) # 创建M折交叉验证的分层抽样
    negfolds <- createFolds(negIND, k = M,list=T) # 创建M折交叉验证的分层抽样
    posfold_values <- lapply(posfolds, function(indices) posIND[indices])
    negfold_values <- lapply(negfolds, function(indices) negIND[indices])
    my_list=c(posfold_values,negfold_values)
    # 获取唯一的名称
    unique_names <- unique(names(my_list))
    # 合并具有相同名称的元素
    merged_list <- lapply(unique_names, function(name) {
      elements <- unname(unlist(my_list[names(my_list) == name]))
      return(elements)
    })
    # 将名称赋给合并后的列表
    names(merged_list) <- unique_names
    folds=merged_list
    #########################################################
    print(length(posfold_values[[1]]))  ##number of positive samples in test set
    print(length(posIND)-length(posfold_values[[1]])) ##number of positive samples in train set
    if(length(posfold_values[[1]])<2){next}
    if((length(posIND)-length(posfold_values[[1]]))<10){next}
    ##########################################################
    data$tag<-as.factor(data$tag)
    featureSet=1:ncol(data)
    
    
    ALLdecision_values=NULL
    ALLprob_values=NULL
    ALLtag=NULL
    sum_SE=0
    sum_SP=0
    sum_ACC=0
    sum_F=0
    sum_MCC=0
    sum_AUC=0
    sum_AUPRC=0
    #循环M次，分别取出第i部分作为测试样本，其余部分作为训练样本
    for(i in 1:M){
      #将数据集分为训练集和测试集
      #i=1
      testIndex <- folds[[i]]
      
      test_data <- data[testIndex, ]
      train_data <- data[-testIndex, ]
  
      test_label=data$tag[testIndex]
      train_label=data$tag[-testIndex]
      
      # test_data<-data[ind==i,featureSet]
      # train_data<-data[ind!=i,featureSet]
      
      #label为样本类别标签，同样选取相应的训练集对应的label(即从所有的label中选取train对应的label
      # test_label=data$tag[ind==i]
      # train_label=data$tag[ind!=i]
      
      
      ###计数测试集里的positive and negative data number###########
      posN=nrow(test_data[test_data$tag==1,])
      negN=nrow(test_data[test_data$tag==0,])
      ########################################################
      ##提取最后一列label以外的feature matrix
      train=train_data[,-ncol(train_data)]   
      test=test_data[,-ncol(test_data)]
      
      ########训练与测试资料标准化######下面的归一化是独立归一化,不合适，应该用train的归一化mode去归一化test
      # maxst <- apply(train, 2, max)
      # minst <- apply(train, 2, min)
      # scaled.train <- as.data.frame(scale(train,center = mins, scale = maxs - mins))
      # maxst <- apply(test, 2, max)
      # minst <- apply(test, 2, min)
      # scaled.test <- as.data.frame(scale(test,center = mins, scale = maxs - mins))
      # scaled.train$tag=train_label
      # scaled.test$tag=test_label
      
      ###############################################################
      
      #数据预处理
      #train_data$tag = as.factor(train_data$tag)
      #test_data$tag = as.factor(test_data$tag)
      
      
      #########normalization###################
      train.min <- apply(train,2, min)
      train.max <- apply(train,2, max)
      
      #########normalization method1正确###################
      train <- as.data.frame(scale(train,center = train.min, scale =train.max - train.min))
      test <- as.data.frame(scale(test,center = train.min, scale =train.max - train.min))

      cols_with_na1 <- colSums(is.na(train)) > 0
      cols_with_na2 <- colSums(is.na(test)) > 0
      cols_with_na <- ifelse(cols_with_na1==T | cols_with_na2==T,TRUE,FALSE)
      #names(train)[cols_with_na]
      train <- train[,!cols_with_na]
      #names(train)[cols_with_na]
      test <- test[,!cols_with_na]
      
      # Lasso筛选变量动态过程图
      #la.eq <- glmnet(x=as.matrix(train), y=as.matrix(train_label), family="gaussian", intercept = F, alpha=1) ##回归
      la.eq <- glmnet(x=as.matrix(train), y=train_label, family="binomial", alpha=1) 
      # 当alpha设置为0则为ridge回归，将alpha设置为0和1之间则为elastic net 
      
      pdf("LASSO_feature_selection.pdf",width=6,height=5)
      #par(mfcol=c(1,1))
      par(cex=1.3)
      #par(mar=c(4,4,2,1))
      #par(mgp = c(2, 0.5, 0))
      plot(la.eq,xvar = "lambda", label = F)
      # 也可以用下面的方法绘制
      #matplot(log(la.eq$lambda), t(la.eq$beta),type="l", main="Lasso", lwd=2)
      dev.off()
      #我们可以看到，当lambda越大，各估计参数相应的也被压缩得更小，而当lambda达到一定值以后，一部分不重要的变量将被压缩为0，代表该变量已被剔除出模型，图中从左至左右断下降的曲线如同被不断增大的lambda一步一步压缩，直到压缩为0。
      
      
      #################用“交叉验证”确定最优lambda###############
      model_cv <- cv.glmnet(x=as.matrix(train), y=train_label, family="binomial", alpha=1, nfolds = 10)
      pdf("LASSO_best_lambda.pdf",width=6,height=5)
      #par(mfcol=c(1,1))
      par(cex=1.3)
      #par(mar=c(4,4,2,1))
      #par(mgp = c(2, 0.5, 0))
      plot(model_cv) 
      dev.off()
      #通过交叉验证，我们可以选择平均误差最小的那个λ，即model_cv$lambda.min,也可以选择平均误差在一个标准差以内的最大的λ，即model_cv$lambda.1se。
      best_lambda <- model_cv$lambda.min
      coef.min <- coef(model_cv,s="lambda.min")
      
      ###########保留至少10个特征###################
      num.selected <- length(rownames(coef.min)[coef.min[, 1] != 0])
      if(num.selected<10){
        ########lambda_values是从大到小排序的#####
        lambda_values <- model_cv$lambda
        num_features <- sapply(lambda_values, function(l) sum(coef(model_cv, s = l) != 0))
        best_lambda <- lambda_values[which.max(num_features >= 10)] #which.max寻找最大值的索引号，这里实际上找到第一满足条件的TRUE
      }
      ##############################################
      #现在可以用最优lambda值构建的最终模型。
      best_model <- glmnet(x=as.matrix(train), y=train_label, family="binomial", alpha=1, lambda = best_lambda)
      coef.best <-coef(best_model)
      
      ###保存模型######
      #saveRDS(best_model,"best_model.rds")
      #best_model <- readRDS("best_model.rds")
      ######################################################################################
      #######################测试集预测###################################################
      y_predicted <- predict(best_model, s = best_lambda,newx = as.matrix(test),decision.values=T)
      #decision.values = TRUE参数会返回决策函数的值，对于二分类问题，这通常是logit变换后的线性预测值，即log(odds ratio)。
      y_predicted2 <- predict(best_model, s = best_lambda,newx = as.matrix(test),type = "response")
      #对于二分类问题，type = "response"参数会返回每个样本属于正类的概率预测。
      y_predicted3 <- predict(best_model, s = best_lambda,newx = as.matrix(test),type = "class")
      #如果你需要将概率转换为类别标签，可以使用type = "class"参数：
      ##########回归预测的评价指标###############
      # #find SST and SSE
      # y=test_label
      # sst <- sum((y - mean(y))^2)
      # sse <- sum((y_predicted - y)^2)
      # 
      # #find R-Squared
      # rsq <- 1 - sse/sst
      ##########################################
      decision_values=y_predicted
      prob_values=y_predicted2
      ###输出预测结果的列联表
      predicted.table=table(test_data$tag,y_predicted3,dnn=c("experiment","prediction"))  # dnn=dimnames names
      ########################################################################################
      ALLdecision_values=c(ALLdecision_values,as.numeric(decision_values))
      ALLprob_values=c(ALLprob_values,as.numeric(prob_values))
      ALLtag=c(ALLtag,as.numeric(as.character(test_data$tag)))
      
      
      ###########################################################################
      # ##为什么prediced.table这里的数据与下面最佳阈值对应的数据不一样,这是因为pre_svm <- predict预测是SVM默认预测，预测用的decision_value的阈值不一定是最佳阈值，预测结果不是最佳的(highest value of SE+SP)
      # TP <-  predicted.table[2,2]
      # FP <- predicted.table[1,2]
      # TN <- predicted.table[1,1]
      # FN <- predicted.table[2,1]
      # SE <- TP / (TP + FN)
      # SP <- TN / (TN + FP)
      # ACC <- (TP + TN) / (TP + TN + FP + FN)
      # TPR<-SE
      # TNR<-SP
      # 
      # alpha=0.5
      # PPV=TP/(TP+FP)   #positive predictive value, also called precision
      # Fmeasure = 1/ (alpha*1/PPV + (1-alpha)*1/TPR)
      # 
      
      
      #classAgreement(predicted.table)
      #confusionMatrix(predicted.table)
      ###################################################################################
      
      
      ######################绘制ROC曲线################################################
      model_roc <- roc(test_data$tag,as.numeric(prob_values))
      # pdf("Lasso_ROC_Prediction1.pdf",width=6,height=5)
      # #par(mfcol=c(1,1))
      # par(cex=1.3)
      # #par(mar=c(4,4,2,1))
      # #par(mgp = c(2, 0.5, 0))
      # plot(model_roc,legacy.axes=T, print.auc=T, auc.polygon=T, max.auc.polygon=T,auc.polygon.col="orange", print.thres=F,main='ROC curve of Lasso')
      # dev.off()
      
      
      scores <- data.frame(predictions=as.numeric(prob_values), labels=as.numeric(as.character(test_data$tag)))
      pr <- pr.curve(scores.class0=scores[scores$labels==1,"predictions"],
                     scores.class1=scores[scores$labels==0,"predictions"],
                     curve=T,max.compute=T, min.compute=T, rand.compute=T)
      
      
      
      ################# 获取最佳阈值,opt里依次包含threshold, specificity, sensitivity###########
      opt <- coords(model_roc, "best",transpose = TRUE)
      # 计算在最佳阈值下混淆矩阵各项的值
      TP <- dim(test_data[(test_data$tag)==1 & prob_values > opt[1],])[1]
      FP <- dim(test_data[(test_data$tag)==0 & prob_values > opt[1],])[1]
      TN <- dim(test_data[(test_data$tag)==0 & prob_values <= opt[1],])[1]
      FN <- dim(test_data[(test_data$tag)==1 & prob_values <= opt[1],])[1]
      
      # 根据混淆矩阵计算特异度、敏感度以及准确度指标
      SE <- TP / (TP + FN) ##也叫召回率Recall
      SP <- TN / (TN + FP)
      ACC <- (TP + TN) / (TP + TN + FP + FN)
      TPR<-SE
      TNR<-SP
      MCC=((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TN+FN)*(TP+FN)*(TN+FP))
      
      alpha=0.5
      PPV=TP/(TP+FP)   #positive predictive value, also called Precision精确率
      Fmeasure = 1/ (alpha*1/PPV + (1-alpha)*1/TPR) #即Fmeasure=2*PPV*SE/(PPV+SE) ##也就是Fmeasure=2*Precision*Recall/(Precision+Recall)
      
      Result=data.frame(SE=SE,SP=SP,ACC=ACC,Fmeasure=Fmeasure,MCC=MCC,AUC=model_roc$auc,AUPRC=pr$auc.integral)
      if(i==1){
        write.table(Result,file=paste0("selected_feature_LASSO_finalResult_averagePrediction_",k,".txt"),quote=FALSE,col.names=T,row.names=F,sep='\t')
      }else{
        write.table(Result,file=paste0("selected_feature_LASSO_finalResult_averagePrediction_",k,".txt"),quote=FALSE,col.names=F,row.names=F,sep='\t',append = T)
      }
      
      sum_SE=sum_SE+SE
      sum_SP=sum_SP+SP
      sum_ACC=sum_ACC+ACC
      sum_F=sum_F+Fmeasure
      sum_MCC=sum_MCC+MCC
      sum_AUC=sum_AUC+model_roc$auc
      sum_AUPRC=sum_AUPRC+pr$auc.integral
      
    }
    
    averSE=sum_SE/M
    averSP=sum_SP/M
    averACC=sum_ACC/M
    averF=sum_F/M
    averMCC=sum_MCC/M
    averAUC=sum_AUC/M
    averAUPRC=sum_AUPRC/M
    
    Result=data.frame(SE=averSE,SP=averSP,ACC=averACC,Fmeasure=averF,MCC=averMCC,AUC=averAUC,AUPRC=averAUPRC)
    write.table(Result,file=paste0("selected_feature_LASSO_finalResult_averagePrediction_",k,".txt"),quote=FALSE,col.names=F,row.names=F,sep='\t',append = T)
    
  }


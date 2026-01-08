rm(list=ls())
gc()
library(data.table)
library(MASS)
library(glmnet)
library(GSEABase)
library(parallel)
library(scalreg)
library(dplyr)
library(tidyr)
library(SIHR)
library(freebird)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('meddic.cpp')
source('meddic.R')

start_time <- Sys.time() 

clinical_data <- fread("TCGA.LUNG_processed_2.csv", header = TRUE)
clinical_data$gender <- factor(clinical_data$gender, levels = c("FEMALE", "MALE"))
clinical_data$gender <- as.numeric(clinical_data$gender)-1

clinical_data$radiation_therapy <- factor(clinical_data$radiation_therapy, levels = c("NO", "YES"))
clinical_data$radiation_therapy <- as.numeric(clinical_data$radiation_therapy)-1

clinical_data$pathologic_stage <- factor(clinical_data$pathologic_stage, levels = c("Stage I","Stage IA","Stage IB","Stage IIA","Stage IIB","Stage IIIA","Stage IIIB","Stage IV"))
clinical_data$pathologic_stage <- as.numeric(clinical_data$pathologic_stage)
clinical_data$pathologic_stage[clinical_data$pathologic_stage<=3] <- 1
clinical_data$pathologic_stage[clinical_data$pathologic_stage>3 & clinical_data$pathologic_stage<=5] <- 2
clinical_data$pathologic_stage[clinical_data$pathologic_stage>5 & clinical_data$pathologic_stage<=7] <- 3
clinical_data$pathologic_stage[clinical_data$pathologic_stage>7] <- 4


# 使用 fread 读取数据 #太强大了不用转csv直接就能用
methylation_data <- fread("HumanMethylation450-2", header = TRUE, sep = "\t")
# 删除包含 NA 的cg
methylation_data_cleaned <- methylation_data[rowSums(is.na(methylation_data)) == 0,]
colnames(methylation_data_cleaned)[colnames(methylation_data_cleaned) == "sample"] <- "cgID"

#把cgID和geneID对应上
mapping_df <- fread("probeMap-illuminaMethyl450_hg19_GPL16304_TCGAlegacy", header = TRUE, sep = "\t")  # 包含 "cgID" 和 "GeneID" 列
colnames(mapping_df)[colnames(mapping_df) == "#id"] <- "cgID"
colnames(mapping_df)[colnames(mapping_df) == "gene"] <- "geneID"
merged_data <- merge(methylation_data_cleaned, mapping_df, by = "cgID")

# 读取 GMT 文件 Molecular Signatures Database
GSEA_data <- getGmt("c5.go.bp.v2024.1.Hs.symbols.gmt.txt")

# 提取gmt基因集，正确使用 GeneSetCollection 对象的方法
GSEA_list <- names(GSEA_data)

# 遍历每个基因集
test_NIE_for_each_geneset <- function(aa){
  p_result <- c()
  set_name <- GSEA_list[aa]
  # 访问当前 GO 基因集
  go_term <- GSEA_data[[set_name]]
  # 提取 geneIds
  GSEA_gene_ids <- geneIds(go_term)
  
  # 拆分GeneID列
  merged_data_separated <- merged_data %>%
    separate_rows(geneID, sep = ",")
  
  # 筛选数据
  select_cpg <- merged_data_separated %>%
    filter(geneID %in% GSEA_gene_ids) %>%
    select(cgID)%>%
    pull(cgID)
  
  # 提取原始数据框中对应的行
  subset_data <- merged_data %>%
    filter(cgID %in% select_cpg)
  
  # design matrix
  df1 <- t(subset_data[,2:908])
  colnames(df1) <- subset_data$cgID
  common_rows <- intersect(rownames(df1), clinical_data$sampleID)
  M <- df1[common_rows,]
  df2 <- clinical_data[clinical_data$sampleID %in% common_rows,]
  X <- as.matrix(df2[,2:5])
  A <- as.matrix(df2[,7])
  Y <- as.matrix(df2[,6])
  d <- ncol(A)
  p <- ncol(M)
  N <- nrow(M)
  X_intercept <- cbind(rep(1,N),X)
  q <- ncol(X_intercept)
  
  if (p <= 50) {
    return(c(set_name, "<50", "<50", "<50", "<50", "<50", "<50"))
  }

    if(p>1000){
      mydat <- data.frame(cbind(Y,A,X,M))
      # 使用SIS进行变量筛选
      # 计算相关系数矩阵
      cor_matrix <- cor(mydat)
      cor_product <- cor_matrix[(d+q+1):(d+q+p),1]*cor_matrix[2,(d+q+1):(d+q+p)]
      top_cpg <- order(abs(cor_product), decreasing = TRUE)[1:1000]
      M <- M[,top_cpg]
      p <- ncol(M)
    }

    Z <- cbind(M,A,X_intercept)
    W <- cbind(A,X_intercept)
    
    vartheta <- (t(M)%*%W)%*%solve(t(W)%*%W)
    hat_gamma <- vartheta%*%rbind(diag(1,d),matrix(0,q,d))
    hat_delta <- solve(t(W)%*%W)%*%t(W)%*%Y
    Sigma_W <- t(W)%*%W/N
    
    #Our proposed method
    dpenalty <- function(lambda, abs_beta){
      if(abs_beta<=lambda){
        return(lambda)
      }else if(abs_beta>lambda&abs_beta<=3.7*lambda){
        return((3.7*lambda-abs_beta)/(3.7-1))
      }else{
        return(0)
      }
    }
    
    #initial estimation-partial lasso
    p.fac <- rep(1,p+d+q)
    p.fac[(p+1):(p+d+q)] <- 0# partial lasso
    #pfit_lasso <- glmnet(Z, Y, alpha = 1, penalty.factor = p.fac, intercept = F)# alpha=1 means lasso
    cv.lasso <- cv.glmnet(Z, Y, alpha = 1, penalty.factor = p.fac, intercept= F)
    theta <- coef(cv.lasso, s= "lambda.1se")
    theta <- matrix(theta)[2:(p+d+q+1)]
    
    #LLA algorithm-HBIC criterion
    LLA.HBIC <- function(lambda){
      
      #update theta
      rr <- 0
      theta1 <- theta
      beta <- theta[1:p]
      p.fac[1:p] <- sapply(abs(beta), function(i) dpenalty(lambda, i))
      pfit_lasso <- glmnet(Z, Y, alpha = 1, penalty.factor = p.fac, intercept = F)# alpha=1 means lasso
      theta <- matrix(coef(pfit_lasso, s = lambda))[2:(p+d+q+1)]
      while(sum((theta-theta1)^2)>0.001&rr<1000){
        rr <- rr+1 #restrict rr to avoid infinite loops
        theta1 <- theta
        beta <- theta[1:p]
        p.fac[1:p] <- sapply(abs(beta), function(i) dpenalty(lambda, i))
        pfit_lasso <- glmnet(Z, Y, alpha = 1, penalty.factor = p.fac, intercept = F)# alpha=1 means lasso
        theta <- matrix(coef(pfit_lasso, s = lambda))[2:(p+d+q+1)]
      }
      
      #HBIC
      df <- sum(theta!=0)
      HBIC <- log(sum((Y-Z%*%theta)^2))+df*log(log(N))*log(p+d+q)/N
      return(list(HBIC,theta,rr))
    }
    
    lambda_set <- seq(0.1,2,by=0.1)
    result <- rep(NA,length(lambda_set))
    #rr_set <- rep(NA,length(lambda_set))
    for (i in 1:length(lambda_set)) {
      result[i] <- LLA.HBIC(lambda_set[i])[[1]]
      #rr_set[i] <- LLA.HBIC(lambda_set[i])[[3]]
    }
    
    num<- which(result==min(result))[1]
    theta.lambda <- LLA.HBIC(lambda_set[num])[[2]]
    beta.lambda <- theta.lambda[1:p]
    alpha1.lambda <- theta.lambda[(p+1):(p+d)]
    alpha2.lambda <- theta.lambda[(p+d+1):(p+d+q)]
    mathcal_A <- which(beta.lambda!=0)
    hat_s <- length(mathcal_A)
    
    #test statistic
    if(hat_s>0 & hat_s<=N){#hat_beta \neq 0
      
      M.A <- M[,mathcal_A]
      Sigma_M.A <- t(M.A)%*%M.A/N
      Sigma_M.A_W <- t(M.A)%*%W/N
      Sigma_W_M.A <- t(W)%*%M.A/N
      
      Tn <- t(hat_gamma)%*%beta.lambda #NIE
      beta.A <- beta.lambda[mathcal_A]
      gamma.A <- hat_gamma[mathcal_A, , drop = F]
      Sigma_A <- t(A)%*%A/N
      Sigma_A_X <- t(A)%*%X_intercept/N
      Sigma_X_A <- t(X_intercept)%*%A/N
      Sigma_X <- t(X_intercept)%*%X_intercept/N
      hat_epsilon1 <- Y-Z%*%theta.lambda
      sigma_epsilon1 <- t(hat_epsilon1)%*%hat_epsilon1/N
      hat_epsilon2 <- M-W%*%t(vartheta)
      hat_epsilon2.A <- hat_epsilon2[,mathcal_A]
      Sigma_epsilon2.A <- t(hat_epsilon2.A)%*%hat_epsilon2.A/N
      if("try-error" %in% class(try(solve(Sigma_M.A-Sigma_M.A_W%*%solve(Sigma_W)%*%Sigma_W_M.A)))){
        p_result[2] <- "singular"
      }else{
        Sigma_T <- c(sigma_epsilon1)*t(gamma.A)%*%solve(Sigma_M.A-Sigma_M.A_W%*%solve(Sigma_W)%*%Sigma_W_M.A)%*%gamma.A/N+solve(Sigma_A-Sigma_A_X%*%solve(Sigma_X)%*%Sigma_X_A)*c(t(beta.A)%*%Sigma_epsilon2.A%*%beta.A)/N
        Vn1 <-t(Tn)%*%solve(Sigma_T)%*%Tn 
        Vn <- Vn1
        p_result[2] <- 1-pchisq(Vn,d)
      }
    }else if(hat_s==0){
      #hat_beta=0
      #if d=1
      if(d==1){
        weight <- function(x){
          ww <- 1/(sqrt(sum(x^2)+1))
          return(ww)
        }
        g <- apply(W,1,weight)
        g <- t(g)
      }else{
        weight <- function(x){
          ww1 <- 1/sqrt(sum(x^2)+1)
          ww2 <- exp(-sum(x^2)/2)
          h <- 1*N^(-0.2)
          ww3 <- prod(dnorm(x/h)/h)#sin(sum(x^2))
          ww <- c(ww1,ww2,ww3)
          return(ww)
        }
        g <- apply(W,1,weight)
      }
      
      Qn <- g%*%(Y-W%*%hat_delta)/sqrt(N)
      Sigma_g <- g%*%t(g)/N
      Sigma_g_W <- g%*%W/N
      Sigma_W_g <- t(W)%*%t(g)/N
      hat_epsilon1 <- Y-W%*%hat_delta
      sigma_epsilon1 <- t(hat_epsilon1)%*%hat_epsilon1/N
      Sigma_Q <- c(sigma_epsilon1)*(Sigma_g-Sigma_g_W%*%solve(Sigma_W)%*%Sigma_W_g)
      Vn2 <- t(Qn)%*%solve(Sigma_Q)%*%Qn
      Vn <- Vn2
      p_result[2] <- 1-pchisq(Vn,d)
    }else{
      p_result[2] <- "singular"
    }
    p_result[1] <- set_name
    p_result[7] <- hat_s
    
    #Guo Xv's method
    if(hat_s>0 & hat_s<=N){
      Z.A <- cbind(A, X_intercept, M.A)
      Sigma_Z.A <- t(Z.A)%*%Z.A/N
      totaleffect <- cbind(diag(1,d),matrix(0,d,q))%*%hat_delta
      NIE.Guo <- totaleffect-alpha1.lambda
      sigma.Guo <- t(Y-W%*%hat_delta)%*%(Y-W%*%hat_delta)/N
      sigma1.Guo <- t(Y-Z%*%theta.lambda)%*%(Y-Z%*%theta.lambda)/(N-hat_s-d-q)
      sigma2.Guo <- max(sigma.Guo-sigma1.Guo,0)
      iden1 <- cbind(diag(d),matrix(0,d,q))
      iden2 <- cbind(diag(d),matrix(0,d,(q+hat_s)))
      #Sigma.Guo <- sigma2.Guo*iden1%*%solve(Sigma_W)%*%t(iden1)+c(sigma1.Guo)*iden1%*%solve(Sigma_W)%*%Sigma_W_M.A%*%solve(Sigma_M.A-Sigma_M.A_W%*%solve(Sigma_W)%*%Sigma_W_M.A)%*%Sigma_M.A_W%*%solve(Sigma_W)%*%t(iden1)
      Sigma.Guo <- sigma2.Guo*iden1%*%solve(Sigma_W)%*%t(iden1)+c(sigma1.Guo)*(iden2%*%solve(Sigma_Z.A)%*%t(iden2)-iden1%*%solve(Sigma_W)%*%t(iden1))
      
       if("try-error" %in% class(try(solve(Sigma.Guo)))){
        p_result[3] <- "singular"
      }else{
        Vn.Guo <- N*t(NIE.Guo)%*%solve(Sigma.Guo)%*%NIE.Guo
        p_result[3] <- 1-pchisq(Vn.Guo,d)
      }
    }else{
      p_result[3] <- "singular"
    }
    
    
    #Lin Yinan's method
    #initial estimator-scalreg
    lasso.Lin <- scalreg(Z, Y, lam0 = "quantile")
    theta.Lin <- coef(lasso.Lin)
    sigma.Lin <- t(Y-W%*%hat_delta)%*%(Y-W%*%hat_delta)/N
    sigma_Z.Lin <- (lasso.Lin$hsigma)^2
    sigma_E.Lin <- max(0,sigma.Lin-sigma_Z.Lin)
    gamma.Lin <- hat_gamma
    beta.Lin <- theta.Lin[1:p]
    NIE.Lin <- t(gamma.Lin)%*%beta.Lin
    #search direction U
    loading.mat <- rbind(gamma.Lin,matrix(0,d+q,d))
    Est <- try(LF(Z, Y, loading.mat, intercept = F, beta.init =  theta.Lin))
    if("try-error" %in% class(Est)){
      p_result[4] <- "singular"
    }else{
      U <- Est$proj.mat
      Z_cen <- scale(Z, center=TRUE, scale=F)
      debiased_NIE.Lin <- NIE.Lin+t(U)%*%t(Z_cen)%*%(Y-Z_cen%*%theta.Lin)/N
      tau <- 1
      li_tran <- cbind(diag(d),matrix(0,d,q))
      Sigma.Lin <- sigma_E.Lin*li_tran%*%solve(Sigma_W)%*%t(li_tran)/N+sigma_Z.Lin*t(U)%*%(t(Z_cen)%*%Z_cen/N)%*%U/N+tau*diag(d)/N
      
      #if d=1
      if(d==1){
        Vn.Lin <- debiased_NIE.Lin/sqrt(Sigma.Lin)
      }else{
        Vn.Lin <- c()
        for (k in 1:d) {
          Vn.Lin[k] <- debiased_NIE.Lin[k]/sqrt(diag(Sigma.Lin)[k])
        }
      }
      p_result[4] <- 2 * (1 - pnorm(max(abs(Vn.Lin)))) * d
      p_result[4] <- min(1, p_result[4])  # less than one after correction
    }
    
    #Zhou Ruixuan's method
    M_star <- M-X_intercept%*%solve(t(X_intercept)%*%X_intercept)%*%t(X_intercept)%*%M
    Y_star <- Y-X_intercept%*%solve(t(X_intercept)%*%X_intercept)%*%t(X_intercept)%*%Y
    
    output.Zhou <- hilma(Y_star, M_star, A, mediation_setting = 'incomplete')
    p_val.Zhou <- output.Zhou$pvalue_beta_hat[1:d,]
    p_result[5] <- p_val.Zhou
    
    #Zhang Qi's method
    
    meddic_out <- meddic_Wrapper(Y,X,A,M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=0,f_score = 'score_sim.rdata')
    inference_out <- meddic_out$results_inference_alpha05
    p_val.Zhang <- inference_out$pval.asymp[,3]
    p_result[6] <- p_val.Zhang
    
  #sfCat(paste("Iteration ", x), sep = "\n")
  return(p_result)
}


cl <- makeCluster(96)

clusterEvalQ(cl, {
  library(data.table)
  library(MASS)
  library(glmnet)
  library(GSEABase)
  library(parallel)
  library(scalreg)
  library(dplyr)
  library(tidyr)
  library(SIHR)
  library(freebird)
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp('meddic.cpp')
  source('meddic.R')
})

clusterExport(cl, c("GSEA_list", 
                    "GSEA_data", 
                    "merged_data", 
                    "clinical_data", 
                    "test_NIE_for_each_geneset"))

result <- parLapply(cl, 1:length(GSEA_list), function(aa) {
  # log_file <- paste0( "geneset_", aa, ".txt")
  # cat(paste("Running iteration", aa, "\n"), file = log_file, append = TRUE)
  test_NIE_for_each_geneset(aa)
})
res_mat <- matrix(unlist(result),byrow= T, ncol = 7)
write.csv(res_mat, file = "2025.9.16_real_example_LUNG.csv")


# #这个是enrich analysis，我们用不到，我们就是retrive from GO
# #ClusterProfiler可以帮助你将基因数据与GO注释关联，从而揭示潜在的生物学机制或差异。
# go_terms <- enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
# go_terms_df <- as.data.frame(go_terms)
# 
# dimm <- c()
# aa <- 0
# for (term in unique(go_terms_df$ID)) {
#   # 获取该GO类别下的基因
#   genes_in_term <- go_terms_df$geneID[go_terms_df$ID == term]
#   
#   # 使用 strsplit 函数分隔基因名称
#   genes_vector <- unlist(strsplit(genes_in_term, "/"))
#   
#   # 从合并数据中筛选相关基因
#   subset_data <- merged_data[merged_data$gene %in% genes_vector, ]
#   aa <- aa+1
#   dimm[aa] <- dim(subset_data)[1]
#   # 保存为CSV文件
#   write.csv(subset_data, paste0("subset_", term, ".csv"))
# }

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

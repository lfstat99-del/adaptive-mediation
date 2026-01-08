rm(list=ls())
gc()
library(parallel)
library(MASS)
library(glmnet)
library(scalreg)
library(SIHR)
library(freebird)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('meddic.cpp')
source('meddic.R')


start_time <- Sys.time() 

testNIE <- function(N, p, q, s, a, c1, c2, d){
  # data generation
  if(a==1){
    A <- mvrnorm(N,rep(0,d),diag(d))
  }else{
    A <- matrix(rbinom(N*d,1,0.5),N,d)
  }
  X <- mvrnorm(N,rep(0,q),diag(q))
  epsilon1 <- rnorm(N,0,1)
  epsilon2 <- mvrnorm(N,rep(0,p),matrix(0.5,p,p)^abs(matrix(1:p,p,p)-t(matrix(1:p,p,p))))
  alpha3 <- matrix(1,p,q)
  gamma <- c1*matrix(rep(c(0.2*(1:s),rnorm(p-s,0,0.1)),d),ncol = d)
  M <- A%*%t(gamma)+X%*%t(alpha3)+epsilon2
  beta <- c2*c(1-0.2*(0:(s-1)),rep(0,p-s))
  alpha1 <- rep(0.5,d)
  alpha2 <- rep(1,q)
  Y <- M%*%beta+A%*%alpha1+X%*%alpha2+epsilon1
  Z <- cbind(M,A,X)
  W <- cbind(A,X)
  
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
  
  lambda_set <- seq(0.05, 0.5, by=0.05)
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
  if(hat_s>0){#hat_beta \neq 0
    
    M.A <- M[,mathcal_A]
    Sigma_M.A <- t(M.A)%*%M.A/N
    Sigma_M.A_W <- t(M.A)%*%W/N
    Sigma_W_M.A <- t(W)%*%M.A/N
    
    Tn <- t(hat_gamma)%*%beta.lambda #NIE
    beta.A <- beta.lambda[mathcal_A]
    gamma.A <- hat_gamma[mathcal_A, , drop = F]
    Sigma_A <- t(A)%*%A/N
    Sigma_A_X <- t(A)%*%X/N
    Sigma_X_A <- t(X)%*%A/N
    Sigma_X <- t(X)%*%X/N
    hat_epsilon1 <- Y-Z%*%theta.lambda
    sigma_epsilon1 <- t(hat_epsilon1)%*%hat_epsilon1/N
    hat_epsilon2 <- M-W%*%t(vartheta)
    hat_epsilon2.A <- hat_epsilon2[,mathcal_A]
    Sigma_epsilon2.A <- t(hat_epsilon2.A)%*%hat_epsilon2.A/N
    Sigma_T <- c(sigma_epsilon1)*t(gamma.A)%*%solve(Sigma_M.A-Sigma_M.A_W%*%solve(Sigma_W)%*%Sigma_W_M.A)%*%gamma.A/N+solve(Sigma_A-Sigma_A_X%*%solve(Sigma_X)%*%Sigma_X_A)*c(t(beta.A)%*%Sigma_epsilon2.A%*%beta.A)/N
    Vn1 <-t(Tn)%*%solve(Sigma_T)%*%Tn 
    Vn <- Vn1
  }else{
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
  }
  res <- ifelse(Vn>qchisq(0.95,d),1,0)
  
  #Guo Xv's method
  if(hat_s>0){
    Z.A <- cbind(A, X, M.A)
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
      res.Guo <- NA
    }else{
      Vn.Guo <- N*t(NIE.Guo)%*%solve(Sigma.Guo)%*%NIE.Guo
      res.Guo <- ifelse(Vn.Guo>qchisq(0.95,d),1,0)
    }
  }else{
    res.Guo <- NA
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
    res.Lin <- NA
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
    res.Lin <- ifelse(max(abs(Vn.Lin))>qnorm(1-0.05/(2*d)),1,0)
  }
  
  
  #Zhou Ruixuan's method
  M_star <- M-X%*%solve(t(X)%*%X)%*%t(X)%*%M
  Y_star <- Y-X%*%solve(t(X)%*%X)%*%t(X)%*%Y
  
  output.Zhou <- hilma(Y_star, M_star, A, mediation_setting = 'incomplete')
  p_val.Zhou <- output.Zhou$pvalue_beta_hat[1:d,]
  res.Zhou <- mean(ifelse(p_val.Zhou<0.05, 1, 0))
  
  
  #Zhang Qi's method
  
  meddic_out <- meddic_Wrapper(Y,X,A,M,thr.lse=0.2,no.penalty.exposure=0,no.penalty.confounder=0,f_score = 'score_sim.rdata')
  inference_out <- meddic_out$results_inference_alpha05
  p_val.Zhang <- inference_out$pval.asymp[,3]
  res.Zhang <- mean(ifelse(p_val.Zhang<0.05, 1, 0))
  
  result_all <- c(res, res.Guo, res.Lin, res.Zhou, res.Zhang)
  
  return(result_all)
}

cl <- makeCluster(96)

clusterEvalQ(cl, {
  library(MASS)
  library(glmnet)
  library(scalreg)
  library(SIHR)
  library(freebird)
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp('meddic.cpp')
  source('meddic.R')
})


q <- 3
p <- 300
N <- 300
s <- 5

for (a in 1:2) {
  for (d in c(1,3)) {
    for (c1 in c(1)) {
      for (c2 in seq(-0.2,0.2,by=0.04)) {
        
        clusterExport(cl, varlist = c('testNIE', 'N', 'p', 'q', 's', 'a', 'c1', 'c2', 'd'))
        
        #clusterExport(cl, varlist = ls())
        
        result <- parLapply(cl, 1:1000, function(i){testNIE(N, p, q, s, a, c1, c2, d)})
        
        res_mat <- matrix(unlist(result),byrow = T,ncol = 5)
        Size1 <- apply(res_mat,2,sum,na.rm =T)/1000
        Size2 <- apply(res_mat,2,mean,na.rm =T)#na.rm =T means remove NA and average the rest
        num_na <- apply(res_mat,2,function(x) sum(is.na(x)))#number of NA
        output <- rbind(Size1,Size2,num_na)
        write.csv(output, file = paste0("med2_H1","_a=",a,"_d=",d,"_c1=",c1,"_c2=",c2,"_n=",N,"_p=",p,".csv"))
      }
    }
  }
}

stopCluster(cl)

end_time <- Sys.time()
end_time - start_time

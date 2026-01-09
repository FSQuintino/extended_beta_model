library(gsl)
library(MASS)
library(latex2exp)
library(AdequacyModel)
library(goftest)

#Carrega os pacotes necessários para realizar o paralelismo
library(foreach)
library(doParallel)
### ----------------------------------------------
# Ajuste dos modelos concorrentes
### ----------------------------------------------

## information criteria
#ll_max: maximum log-likelihood
#n_par: number of parameters to be estimated
#n: sample size
aic<-function(ll_max, n_par){-2*ll_max+2*n_par}
bic<-function(ll_max, n_par,n){-2*ll_max+n_par*log(n) }
edc<-function(ll_max, n_par,n){-2*ll_max+log(log(n))*n_par}




# Beta Distribution: PDF and CDF
beta_pdf <- function(par, x) {
  alpha=par[1]
  beta=par[2]
  dbeta(x, shape1 = alpha, shape2 = beta)
}

beta_cdf <- function(par, x) {
  alpha=par[1]
  beta=par[2]
  pbeta(x, shape1 = alpha, shape2 = beta)
}




# Kumaraswamy Distribution: PDF and CDF
pdf_kw <- function(par,x){
  a = par[1]
  b = par[2]
  (a*b*(x^(a-1))*(1-x^a)^(b-1))
}
cdf_kw <- function(par,x){
  a = par[1]
  b = par[2]
  1 - (1 - x^a)^b
}


# bimodal Beta - Probability density function.
#?beta#OK
# incomplete beta function ratio
library(zipfR)
#?Ibeta#incomplete beta
#?Rbeta#incomplete beta regularized 
pdf_bbeta <- function(par,x){
  a = par[1]
  b = par[2]
  r=par[3]
  d=par[4]
  z=1+r-2*d*(a/(a+b))+(d^2)*(a*(a+1))/((a+b)*(a+b+1))
  c1=(r+(1-d*x)^2)/(z*beta(par[1],par[2]))
  return(ifelse(x>0 & x<=1, c1*(x^(a-1))*(1-x)^(b-1), 0))
}
cdf_bbeta <- function(par,x){
  a = par[1]
  b = par[2]
  r=par[3]
  d=par[4]
  
  z=1+r-2*d*(a/(a+b))+(d^2)*(a*(a+1))/((a+b)*(a+b+1))
  
  c1=(1+r)*Rbeta(x,a,b) - (2*d*Ibeta(x,a,b)/beta(a,b)) + (d^2)*Ibeta(x,a+2,b)/beta(a,b)
  
  return(  ifelse(x<0, 0,
                  ifelse(x>1, 1,c1/z)))
}
#cdf_bbeta(c(1,1,0.5,1), runif(100))


### ---------------------------------------------
### Ajustes dos modelos beta e kumaraswamy
### ---------------------------------------------
#library(AdequacyModel)
fit_beta=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(20,20)){
  result_1 = goodness.fit(pdf = beta_pdf, cdf = beta_cdf,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}

fit_kumaraswamy=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(10,10)){
  result_1 = goodness.fit(pdf = pdf_kw, cdf = cdf_kw,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}




# Data reported in Dasgupta (2011):
d10=c(0.04, 0.02, 0.06, 0.12, 0.14, 0.08, 0.22, 0.12, 0.08, 0.26, 0.24, 0.04, 0.14, 0.16,
      0.08, 0.26, 0.32, 0.28, 0.14, 0.16, 0.24, 0.22, 0.12, 0.18, 0.24, 0.32, 0.16, 0.14,
      0.08, 0.16, 0.24, 0.16, 0.32, 0.18, 0.24, 0.22, 0.16, 0.12, 0.24, 0.06, 0.02, 0.18,
      0.22, 0.14, 0.06, 0.04, 0.14, 0.26, 0.18, 0.16)


set.seed(1)
p10b=fit_beta(d10)
p10k=fit_kumaraswamy(d10)
p10ebb=ebb_fit(d10, p10b$par)




library(e1071)
(tab=(D1=data.frame(
  Dataset=c("W"),
  n=c( length(d10)),
  Min=c(min(d10)),
  fStQu=c( quantile(d10,0.25)),
  Median=c(median(d10)),
  Mean=c( mean(d10)),
  trdQu=c( quantile(d10,0.75)),
  Max=c( max(d10)), 
  Sd=c(sd(d10)),
  CS=c( skewness(d10)),
  CK=c(kurtosis(d10))     
)))

write.csv2(tab, "summary_datasets.csv", row.names = F)

library("xtable")
xtable(tab)

###
dir()

readLines("rfam08.csv",10)
readLines("risfam08.csv" ,10)

readLines("rfam08.csv",10)
readLines("risfam08.csv" ,10)

rfam08=read.csv2("rfam08.csv", sep = ",")
risfam08=read.csv2("risfam08.csv", sep = ",")

is.numeric(rfam08$Y)
rfam08$Y=as.numeric(rfam08$Y)
is.numeric(risfam08$C)
risfam08$C=as.numeric(risfam08$C)

rfam08$nquest
risfam08$nquest
names(rfam08)%in%names(risfam08)
names(risfam08)%in%names(rfam08)
names(rfam08)[names(rfam08)=="Y"]="Y_fam"
dfam=merge.data.frame(rfam08,risfam08, by="nquest")

#verificando se ha valores negativos
#despesa
sum(dfam$C<0)
#receita
sum(dfam$Y_fam<0)

#Removendo valores negativos
dfam2=dfam[dfam$C>0 & dfam$Y_fam>0, names(dfam)%in%c("C", "Y_fam")]
sum(dfam2$Y_fam>dfam2$C)/nrow(dfam2)

cor(dfam2$Y_fam, dfam2$C)
#0.7165915

d9=dfam2$Y_fam/(dfam2$Y_fam + dfam2$C)
# hist(d9)
# hist(d9, breaks = 50)
# hist(d9, breaks = 50, probability = T)

set.seed(1)
p9b=fit_beta(d9, lim_sup = c(30,30))
p9k=fit_kumaraswamy(d9, lim_sup = c(30,30))
#p9ebb=ebb_fit(d9, p9b$par) MUITO LENTO

p9ebb = estimador_rho_gamma3(x=dfam2$Y_fam/1000, y=dfam2$C/1000,
                             alpha1=p9b$par[1], alpha2=p9b$par[2])
#Fiz tentativas e erro com 10^(-k) para escolher o denominador da escala



## Histogramas e ECDFs


## d9
x11()
par(mfrow=c(1,2))

#x_points=seq(0,1,by=1/100)

pdf_eb9 = deb_fast(x_points, par=c(p9b$par, p9ebb$rho_hat))

hist(d9, probability = TRUE, breaks = 15,
     main="", 
     #main=TeX("$Z_1$"),
     xlab = "z", 
     cex.main = 1.5)#aumentar o tamanho do titulo)
lines(x_points, beta_pdf(x_points, par=p9b$par), col="green")
lines(x_points, pdf_kw(x_points, par=p9k$par), col="blue")
lines(x_points, pdf_eb9, col="red")
#-0.1648569

legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("Beta", "Kumaraswamy", "EBB"), 
       lty=c(1, 1, 1), #pch=15:16,
       col=c("green", "blue","red")) 


#x_points2=seq(0,1, length.out=30)
cdf9b=beta_cdf(x_points2, par=p9b$par)
cdf9k=cdf_kw(x_points2, par=p9k$par)
cdf9eb=peb_integral(x_points2, par=c(p9b$par, p9ebb$rho_hat))


plot(ecdf(d9), main="",
     #main=TeX("$Z_1$"),
     xlab = "z", ylab = "Fn(z)")
lines(x_points2, cdf9b, col="green")
lines(x_points2, cdf9k, col="blue")
lines(x_points2, cdf9eb, col="red")

legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF", "Beta", "Kumaraswamy", "EBB"), 
       lty=c(1,1, 1,1), #pch=15:16,
       col=c("black","green", "blue","red")) 




par_est_z = data.frame(RV=rep("Z",3),
                       Model=c("Beta", "Kumaraswamy", "EBB"),
                       Parameter1=c(p9b$par[1], p9k$par[1], p9b$par[1]),
                       Parameter2=c(p9b$par[2], p9k$par[2], p9b$par[2]),
                       Parameter3=c(NA, NA, p9ebb$rho_hat)
)


llf_d9=sum(log(deb_fast(d9, par=c(p9b$par, p9ebb$rho_hat))))
(results_z=data.frame(
  #Model=c("Beta", "Kumaraswamy", "EBB"),
  llf=c(p9b$llf_max, p9k$llf_max, llf_d9),
  AIC=c(p9b$AIC, p9k$AIC, aic(llf_d9, 3)),
  BIC=c(p9b$BIC, p9k$BIC, bic(llf_d9, 3, length(d9)))
))

(r9=cbind.data.frame(par_est_z, results_z))



## d10
x11()
par(mfrow=c(1,2))

x_points=seq(0,1,by=1/100)

pdf_eb10 = deb_fast(x_points, par=p10ebb$par)

hist(d10, probability = TRUE, breaks = 5,
     main="", ylim=c(0,5),
     #main=TeX("$Z_1$"),
     xlab = "z", 
     cex.main = 1.5)#aumentar o tamanho do titulo)
lines(x_points, beta_pdf(x_points, par=p10b$par), col="green")
lines(x_points, pdf_kw(x_points, par=p10k$par), col="blue")
lines(x_points, pdf_eb10, col="red")

legend("topright", inset=.1,
       #y=.3,
       #title="",
       c("Beta", "Kumaraswamy", "EBB"), 
       lty=c(1, 1, 1), #pch=15:16,
       col=c("green", "blue","red")) 



lines(x_points, beta_pdf(x_points, par=c(2.80982, 11.84246)), col="magenta")
lines(x_points,deb_fast(x_points, par=c(2.8098200, 11.8424600,  0.6782148)), col="purple")
#ebb_fit(d10, c(2.80982, 11.84246))


x_points2=seq(0,1, length.out=20)
cdf10b=beta_cdf(x_points2, par=p10b$par)
cdf10k=cdf_kw(x_points2, par=p10k$par)
cdf10eb=peb_integral(x_points2, par=p10ebb$par)

cdf10b2=beta_cdf(x_points2, par=c(2.8098200, 11.8424600))
cdf10eb2=peb_integral(x_points2, par=c(2.8098200, 11.8424600,  0.6782148))


plot(ecdf(d10), main="",
     #main=TeX("$Z_1$"),
     xlab = "z", ylab = "Fn(z)")
lines(x_points2, cdf10b, col="green")
lines(x_points2, cdf10k, col="blue")
lines(x_points2, cdf10eb, col="red")

lines(x_points2, cdf10b2, col="magenta")
lines(x_points2, cdf10eb2, col="purple")

estim_beta_ecdf(x)



legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF", "Beta", "Kumaraswamy", "EBB"), 
       lty=c(1,1, 1,1), #pch=15:16,
       col=c("black","green", "blue","red")) 



par_est_w = data.frame(RV=rep("W",3),
                       Model=c("Beta", "Kumaraswamy", "EBB"),
                       Parameter1=c(p10b$par[1], p10k$par[1], p10ebb$par[1]),
                       Parameter2=c(p10b$par[2], p10k$par[2], p10ebb$par[2]),
                       Parameter3=c(NA, NA, p10ebb$par[3])
)

(results_w=data.frame(
  #Model=c("Beta", "Kumaraswamy", "EBB"),
  llf=c(p7b$llf_max, p7k$llf_max, p7ebb$llf_max),
  AIC=c(p7b$AIC, p7k$AIC, p7ebb$AIC),
  BIC=c(p7b$BIC, p7k$BIC, p7ebb$BIC)
))

(r10=cbind.data.frame(par_est_w, results_w))


(resultados=rbind.data.frame(r1, r7, r9, r10))

xtable(resultados)


###-------------------------------------------
### Comparacao de metodos de estimacao
###-------------------------------------------


p10ebb$par
p10b_kernel=estim_beta_kernel(d10)
p10b_ecdf=estim_beta_ecdf(d10)
p10b_map=beta_bayes_map(d10)

(p10ebb_kernel=ebb_fit(d10, c(p10b_kernel$alpha, p10b_kernel$beta)))
(p10ebb_ecdf=ebb_fit(d10, c(p10b_ecdf$alpha, p10b_ecdf$beta)))
(p10ebb_map=ebb_fit(d10, c(p10b_map$alpha, p10b_map$beta)))


# resultados_comparison = data.frame(
#   Method=c("MLE", "Kernel", "ECDF", "Bayes MAP"),
#   Alpha=c(p10b$par[1], p10b_kernel$alpha, p10b_ecdf$alpha, p10b_map$alpha),
#   Beta=c(p10b$par[2], p10b_kernel$beta, p10b_ecdf$beta, p10b_map$beta),
#   Rho_hat=c(p10ebb$par[3], p10ebb_kernel$par[3], p10ebb_ecdf$par[3], p10ebb_map$par[3]),
#   LogLik=c(p10b$llf_max,
#            sum(log(deb_fast(d10, par=c(p10b_kernel$alpha, p10b_kernel$beta, p10ebb_kernel$par[3])))),
#            sum(log(deb_fast(d10, par=c(p10b_ecdf$alpha, p10b_ecdf$beta, p10ebb_ecdf$par[3])))),
#            sum(log(deb_fast(d10, par=c(p10b_map$alpha, p10b_map$beta, p10ebb_map$par[3]))))
#            )
# )
resultados_comparison = data.frame(
  Method=c("MLE", "Kernel", "ECDF", "Bayes MAP"),
  Alpha=c(p10b$par[1], p10b_kernel$alpha, p10b_ecdf$alpha, p10b_map$alpha),
  Beta=c(p10b$par[2], p10b_kernel$beta, p10b_ecdf$beta, p10b_map$beta),
  Rho_hat=c(p10ebb$par[3], p10ebb_kernel$par[3], p10ebb_ecdf$par[3], p10ebb_map$par[3]),
  LogLik=c(sum(log(deb_fast(d10, par=c(p10b$par[1], p10b$par[2], p10ebb$par[3])))),
           sum(log(deb_fast(d10, par=c(p10b_kernel$alpha, p10b_kernel$beta, p10ebb_kernel$par[3])))),
           sum(log(deb_fast(d10, par=c(p10b_ecdf$alpha, p10b_ecdf$beta, p10ebb_ecdf$par[3])))),
           sum(log(deb_fast(d10, par=c(p10b_map$alpha, p10b_map$beta, p10ebb_map$par[3]))))
  )
)


#kolmogorov-Smirnov test
#?ks.test
ks_result_d10 = ks.test(d10, FZ4, par=p10ebb$par,
                        simulate.p.value = TRUE, B = 1000)
ks_result_d10$statistic
ks_result_d10$p.value

ks_result_d10_kernel = ks.test(d10, FZ4, par=p10ebb_kernel$par,
                               simulate.p.value = TRUE, B = 1000)

ks_result_d10_ecdf = ks.test(d10, FZ4, par= p10ebb_ecdf$par,
                             simulate.p.value = TRUE, B = 1000)

ks_result_d10_map = ks.test(d10, FZ4, par=p10ebb_map$par,
                            simulate.p.value = TRUE, B = 1000)


## teste Anderson-Darling
#library(goftest)
?ad.test
ad_result_d10 = ad.test(d10, FZ4, par=p10ebb$par,
                        simulate.p.value = TRUE,
                        B = 1000)
ad_result_d10$statistic
ad_result_d10$p.value


ad_result_d10_kernel = ad.test(d10, FZ4, par=p10ebb_kernel$par,
                               simulate.p.value = TRUE,
                               B = 1000)

ad_result_d10_ecdf = ad.test(d10, FZ4, p10ebb_ecdf$par,
                             simulate.p.value = TRUE,
                             B = 1000)

ad_result_d10_map = ad.test(d10, FZ4, par=p10ebb_map$par,
                            simulate.p.value = TRUE,
                            B = 1000)

#?cvm.test
cvm_result_d10 = cvm.test(d10, FZ4, par=p10ebb$par,
                          simulate.p.value = TRUE,
                          B = 1000)
cvm_result_d10$statistic
cvm_result_d10$p.value


cvm_result_d10_kernel = cvm.test(d10, FZ4, par=p10ebb_kernel$par,
                                 simulate.p.value = TRUE,
                                 B = 1000)

cvm_result_d10_ecdf = cvm.test(d10, FZ4, p10ebb_ecdf$par,
                               simulate.p.value = TRUE,
                               B = 1000)

cvm_result_d10_map = cvm.test(d10, FZ4, par=p10ebb_map$par,
                              simulate.p.value = TRUE,
                              B = 1000)


(resultados_tests = data.frame(
  Method=c("MLE", "Kernel", "ECDF", "Bayes MAP"),
  KS_Statistic=c(ks_result_d10$statistic,
                 ks_result_d10_kernel$statistic,
                 ks_result_d10_ecdf$statistic,
                 ks_result_d10_map$statistic),
  KS_p_value=c(ks_result_d10$p.value,
               ks_result_d10_kernel$p.value,
               ks_result_d10_ecdf$p.value,
               ks_result_d10_map$p.value),
  AD_Statistic=c(ad_result_d10$statistic,
                 ad_result_d10_kernel$statistic,
                 ad_result_d10_ecdf$statistic,
                 ad_result_d10_map$statistic),
  AD_p_value=c(ad_result_d10$p.value,
               ad_result_d10_kernel$p.value,
               ad_result_d10_ecdf$p.value,
               ad_result_d10_map$p.value),
  CvM_Statistic=c(cvm_result_d10$statistic,
                  cvm_result_d10_kernel$statistic,
                  cvm_result_d10_ecdf$statistic,
                  cvm_result_d10_map$statistic),
  CvM_p_value=c(cvm_result_d10$p.value,
                cvm_result_d10_kernel$p.value,
                cvm_result_d10_ecdf$p.value,
                cvm_result_d10_map$p.value)
))


resultados_comparacao_goftest=cbind.data.frame(resultados_comparison, resultados_tests)
write.csv2(resultados_comparacao_goftest,
           "resultados_comparacao_metodos_ebbeta_d10.csv",
           row.names = F)

x11()
par(mfrow=c(1,2))

x_points=seq(0,1,by=1/100)

#pdf_eb10 = deb_fast(x_points, par=p10ebb$par)
pdf_eb10_kernel = deb_fast(x_points, par=p10ebb_kernel$par)
pdf_eb10_ecdf = deb_fast(x_points, par=p10ebb_ecdf$par)
pdf_eb10_map = deb_fast(x_points, par=p10ebb_map$par)


hist(d10, probability = TRUE, breaks = 5,
     main="", ylim=c(0,6),
     #main=TeX("$Z_1$"),
     xlab = "z", 
     cex.main = 1.5)#aumentar o tamanho do titulo)

lines(x_points, pdf_eb10, col="red", lty=11)
lines(x_points, pdf_eb10_kernel, col="magenta", lty=2)
lines(x_points,pdf_eb10_ecdf, col="purple", lty=3)
lines(x_points,pdf_eb10_map, col="brown", lty=4)




legend("topright", inset=.1,
       #y=.3,
       title="Method",
       c("MLE", "Kernel", "ECDF-based", "MAP"), 
       lty=c(1, 2, 3, 4), #pch=15:16,
       col=c("red", "magenta","purple", "brown")) 





x_points2=seq(0,1, length.out=30)
#cdf10eb=peb_integral(x_points2, par=p10ebb$par)
cdf10eb_kernel=peb_integral(x_points2, par=p10ebb_kernel$par)
cdf10eb_ecdf=peb_integral(x_points2, par=p10ebb_ecdf$par)
cdf10eb_map=peb_integral(x_points2, par=p10ebb_map$par)

plot(ecdf(d10), main="",
     #main=TeX("$Z_1$"),
     xlab = "z", ylab = "Fn(z)")

lines(seq(0,1, length.out=20), cdf10eb, col="red", lty=1)
lines(x_points2, cdf10eb_kernel, col="magenta", lty=2)
lines(x_points2, cdf10eb_ecdf, col="purple", lty=3)
lines(x_points2, cdf10eb_map, col="brown", lty=4)


legend("topleft", inset=.1,
       #y=.3,
       title="Method",
       c("ECDF", "MLE", "Kernel", "ECDF-based", "MAP"), 
       lty=c(1,1, 2,3,4), #pch=15:16,
       col=c("black","red", "magenta","purple", "brown")) 





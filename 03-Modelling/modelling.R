library(gsl)
library(MASS)
library(latex2exp)
library(AdequacyModel)
library(goftest)

#Carrega os pacotes necessários para realizar o paralelismo
library(foreach)
library(doParallel)

### -------------------------------------------
### reading auxiliar functions
### -------------------------------------------

setwd("G:/02 - Auxiliar functions")
source("auxiliar_functions.R")

### -------------------------------------------
### reading dataset
### -------------------------------------------

setwd("G:/01 - Dataset")

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
# p9k=fit_kumaraswamy(d9, lim_sup = c(30,30))

#calcular o tempo computacional
#tempo inicial
start_time <- Sys.time()
set.seed(1)
p9b=fit_beta(d9, lim_sup = c(30,30))
p9ebb = estimador_rho_gamma3(x=dfam2$Y_fam/1000, y=dfam2$C/1000,
                             alpha1=p9b$par[1], alpha2=p9b$par[2])
#Fiz tentativas e erro com 10^(-k) para escolher o denominador da escala
#tempo final
end_time <- Sys.time()
#tempo computacional
(ebb_pmle = end_time - start_time)
#Time difference of 50.80358 secs


### ----------------------------------------------------
### Comparison of two-steps (pseudolikelihood) and MLE
### ----------------------------------------------------
(param_est = c(p9b$par, p9ebb$rho_hat))
# [1] 16.09642 12.84415  0.78111

(llf_pmle=sum(log(deb_fast(x=d9, param_est))))
# [1] 7922.921
#llf_pmle=sum(log(f_Z(x=d9, param_est)))

aic(llf_pmle, 3)
bic(llf_pmle, 3, length(d9))

llf = function(data, param){
  -sum(log(deb_fast(data, param)))
}

#Otimizar llf
set.seed(1)
start_time_mle <- Sys.time()
opt_llf = optim(par=param_est,
                fn=function(par) -llf(d9, par),
                method="L-BFGS-B",
                lower=c(0.01, 0.01, -1), upper=c(20, 20, 1))
#tempo final
end_time_mle <- Sys.time()
#tempo computacional
(ebb_mle = end_time_mle - start_time_mle)



ebb_pmle; ebb_mle
#converter ambos para minutos
ebb_pmle_minutes = as.numeric(ebb_pmle, units = "mins")
ebb_mle_minutes = as.numeric(ebb_mle, units = "mins")

param_est
opt_llf$par
#[1] 0.01 20.00 1.00 #parametros batendo nas extremidades da caixa =(

opt_llf$value
#[1] -195973.8 

# opt_llf$convergence
# [1] 0
### Funcao aprimorou muito o log-verossimilhança,
## mas os parâmetros estão nas extremidades da caixa de otimização,
# o que pode indicar um problema de convergência ou
# que os limites impostos são muito restritivos.
## Seria interessante revisar os limites ou a função de otimização para garantir
# que os parâmetros possam ser estimados adequadamente.



#Otimizar llf
set.seed(1)
start_time_mle2 <- Sys.time()
opt_llf = optim(par=param_est,
                fn=function(par) -llf(d9, par),
                method="L-BFGS-B",
                lower=c(10^(-4), 10^(-4), -1-10^(-4)), upper=c(40, 40, 1-10^(-4)))
#tempo final
end_time_mle2 <- Sys.time()
#tempo computacional
(ebb_mle2 = end_time_mle2 - start_time_mle2)



#Fiz tentativas e erro com 10^(-k) para escolher o denominador da escala
#tempo final
set.seed(1)
p9k=fit_kw(d9, lim_sup = c(30,30))


# set.seed(1)
# p9tn=fit_tnormal(d9) MUITO LENTO
set.seed(1)
p9ln=fit_logitnormal(d9) 
set.seed(1)
p9tl=fit_toppleone(d9)
set.seed(1)
p9ul=fit_ulindley(d9)
set.seed(1)
p9ug=fit_ugompertz(d9)
set.seed(1)
p9uig=fit_uinvgauss(d9)
set.seed(1)
p9bb=fit_bbeta(d9)
set.seed(1)
p9uw=fit_uweibull(d9)
set.seed(1)
p9ubbs=fit_ubbs(d9)
set.seed(1)
p9uf=fit_uf(d9)




#montar tabela com resultados
tabela_resultados=data.frame(Modelo=c("Beta", "Kumaraswamy", "Logit-Normal", #"Truncated Normal",
                                      "Topp-Leone", "Unit-Lindley", "Unit-Gompertz", "Unit-Inverse Gaussian",
                                      "Bimodal Beta", "Unit-Weibull", "UBBS", "Unit-Frechet"),
                             Estimates =c(paste(round(p9b$par[1], 4),round(p9b$par[2], 4), sep = ", "),
                                          paste(round(p9k$par[1], 4),round(p9k$par[2], 4), sep = ", "),
                                          #paste(round(p9tn$par[1], 4),round(p9tn$par[2], 4), sep = ", "),
                                          paste(round(p9ln$par[1], 4),round(p9ln$par[2], 4), sep = ", "),
                                          paste(round(p9tl$par[1], 4), sep = ", "),
                                          paste(round(p9ul$par[1], 4), sep = ", "),
                                          paste(round(p9ug$par[1], 4),round(p9ug$par[2], 4), sep = ", "),
                                          paste(round(p9uig$par[1], 4),round(p9uig$par[2], 4), sep = ", "),
                                          paste(round(p9bb$par[1], 4),round(p9bb$par[2], 4),
                                                round(p9bb$par[3], 4), round(p9bb$par[4], 4),
                                                sep = ", "),
                                          paste(round(p9uw$par[1], 4),round(p9uw$par[2], 4), sep = ", "),
                                          paste(round(p9ubbs$par[1], 4),round(p9ubbs$par[2], 4),
                                                round(p9ubbs$par[3], 4),
                                                sep = ", "),
                                          paste(round(p9uf$par[1], 4),round(p9uf$par[2], 4),
                                                round(p9uf$par[3], 4),
                                                sep = ", ")),
                             LLF=c(p9b$llf_max, p9k$llf_max, #p9tn$llf_max,
                                   p9ln$llf_max,
                                   p9tl$llf_max, p9ul$llf_max, p9ug$llf_max, p9uig$llf_max,
                                   p9bb$llf_max, p9uw$llf_max, p9ubbs$llf_max, p9uf$llf_max),
                             AIC=c(p9b$AIC, p9k$AIC, #p9tn$AIC, 
                                   p9ln$AIC,
                                   p9tl$AIC, p9ul$AIC, p9ug$AIC, p9uig$AIC,
                                   p9bb$AIC, p9uw$AIC, p9ubbs$AIC, p9uf$AIC),
                             BIC=c(p9b$BIC, p9k$BIC, #p9tn$BIC,
                                   p9ln$BIC,
                                   p9tl$BIC, p9ul$BIC, p9ug$BIC, p9uig$BIC,
                                   p9bb$BIC, p9uw$BIC, p9ubbs$BIC, p9uf$BIC))

setwd("G:/Meu Drive/Pesquisa/Submetidas/EBeta/Round 1/R/comparações")
write.csv2(tabela_resultados, "tabela_resultados.csv", row.names = F)
library(xtable)
print(xtable(tabela_resultados), include.rownames = F)

### ----------------------------------------
### Analises graficas
### ----------------------------------------
#fit densities
x_points=seq(0,1,length.out=100)[-c(1,100)]
pdf_eb9 = deb_fast(x_points, par=c(p9b$par, p9ebb$rho_hat))

db = beta_pdf(x_points, par=p9b$par)
dkw = pdf_kw(x_points, par=p9k$par)
dtl = pdf_toppleone(x_points, par=p9tl$par)
dul=pdf_ulindley(x_points, par=p9ul$par)
dln = pdf_logitnormal(x_points, par=p9ln$par)
dug = pdf_ugompertz(x_points, par=p9ug$par)
duig = pdf_uinvgauss(x_points, par=p9uig$par)
duw = pdf_uweibull(x_points, par=p9uw$par)
duf = pdf_uf(x_points, par=p9uf$par)
dubbs = pdf_ubbs(x_points, par=p9ubbs$par)
dbb = pdf_bbeta(x_points, par=p9bb$par)

x11()
hist(d9, probability = TRUE, breaks = 15,
     main="", 
     #main=TeX("$Z_1$"),
     ylim=c(0, 8),
     xlab = "z", 
     cex.main = 1.5)#aumentar o tamanho do titulo)


#1 parameter
lines(x_points,dtl,col="black", lty=1)
lines(x_points, dul, col="black", lty=2)
#2 parameters
lines(x_points, db, col="blue", lty=1)
lines(x_points, dkw, col="blue", lty=2)
lines(x_points, dln, col="blue", lty=3)
#lines(x_points, pdf_tnormal(x_points, par=p9tn$par), col="purple")
lines(x_points, dug, col="blue", lty=4)
lines(x_points, duig, col="blue", lty=5)
lines(x_points, duw, col="blue", lty=6)
# 3 parameters
lines(x_points, duf, col="magenta", lty=1)
lines(x_points, dubbs, col="magenta", lty=2)
# 4 parameters
lines(x_points,dbb, col="green")
# EBB model
lines(x_points, pdf_eb9, col="red", lwd=2)

legend("topleft", #inset=.1,
       #y=.3,
       #title="",
       c("EB", "Topp-Leone", "Unit-Lindley",
         "Beta",  
         "Kumaraswamy", "Logit-Normal", "Unit-Gompertz",  "Unit-Inverse Gaussian", "Unit-Weibull", "Unit-Frechet",
         "UBBS", "Bimodal Beta"), 
       lty=c(1,1,2,1,2,3,4,5,6,1,2,1), #pch=15:16,
       col=c("red","black","black","blue","blue","blue","blue","blue","blue",
             "magenta","magenta","green"), 
       cex = 1.5)


x_points2=seq(0,1, length.out=30)
cdf9b=beta_cdf(x_points2, par=p9b$par)
cdf9k=cdf_kw(x_points2, par=p9k$par)
cdf9eb=peb_integral(x_points2, par=c(p9b$par, p9ebb$rho_hat))
cdf9ln=cdf_logitnormal(x_points2, par=p9ln$par)
#cdf9tn=cdf_tnormal(x_points2, par=p9tn$par)
cdf9tl=cdf_toppleone(x_points2, par=p9tl$par)
cdf9ul=cdf_ulindley(x_points2, par=p9ul$par)
cdf9ug=cdf_ugompertz(x_points2, par=p9ug$par)
cdf9uig=cdf_uinvgauss(x_points2, par=p9uig$par)
cdf9bb=cdf_bbeta(x_points2, par=p9bb$par)
cdf9uw=cdf_uweibull(x_points2, par=p9uw$par)
cdf9ubbs=cdf_ubbs(x=x_points2, par=p9ubbs$par)
#cdf9ubbs[29]=1
cdf9uf=cdf_uf(x_points2, par=p9uf$par)


x11()
par(mfrow=c(4,3))
plot(ecdf(d9),
     main="EB", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9eb, col="blue", lwd=2)
plot(ecdf(d9), 
     main="Topp-Leone", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9tl, col="blue", lwd=2)

plot(ecdf(d9),
     main="Unit-Lindley", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9ul, col="blue", lwd=2)

plot(ecdf(d9), 
     main="Beta", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9b, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 

plot(ecdf(d9), 
     main="Kumaraswammy", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9k, col="blue", lwd=2)




plot(ecdf(d9), 
     main="Logit-Normal", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9ln, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 





plot(ecdf(d9), 
     main="Unit-Gompertz", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9ug, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 

plot(ecdf(d9), 
     main="Unit-Inverse Gaussian", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9uig, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 


plot(ecdf(d9), 
     main="Unit-Weibull", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9uw, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 



plot(ecdf(d9),
     main="Unit-Frechet", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9uf, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 

plot(ecdf(d9), 
     main="UBBS", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9ubbs, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 

plot(ecdf(d9),
     main="Bimodal Beta", cex.main = 2,
     xlab = "z", ylab = "F(z)", lwd=2)
lines(x_points2, cdf9bb, col="blue", lwd=2)
legend("topleft", inset=.1,
       #y=.3,
       #title="",
       c("ECDF"), 
       lty=c(1), #pch=15:16,
       col=c("black")) 



#### ---------------------------------------------
### Normal QQPlot
### ---------------------------------------------


# x11()
# par(mfrow=c(3,4))

dados=sort(d9)
# Empirical CDF
ecdf_func <- ecdf(dados)  # Empirical CDF function
empirical_cdf <- ecdf_func(dados)  # Values of the empirical CDF at data points


#### EBB
#rebb=qnorm(peb_integral(dados, par=c(p9b$par, p9ebb$rho_hat)))
rebb_finite=rebb[is.finite(rebb)]
# qqnorm(rebb_finite, main = "EBB")
# qqline(rebb_finite, col = "red", lwd = 2)
#write.csv2(rebb_finite, file="rebb.csv")


#### "Topp-Leone"
rtl=qnorm(cdf_toppleone(dados, par=p9tl$par))
# qqnorm(rtl, main = "Topp-Leone")
# qqline(rtl, col = "red", lwd = 2)

#### "Unit-Lindley"
rul=qnorm(cdf_toppleone(dados, par=p9ul$par))
# qqnorm(rul, main = "Unit-Lindley")
# qqline(rul, col = "red", lwd = 2)

#### Beta
rb=qnorm(beta_cdf(dados, par=p9b$par))
# qqnorm(rb, main = "Beta")
# qqline(rb, col = "red", lwd = 2)


#### Kumaraswamy
rkw=qnorm(cdf_kw(dados, par=p9k$par))
# qqnorm(rkw, main = "Kumaraswamy")
# qqline(rkw, col = "red", lwd = 2)


#### "Logit-Normal"
rln=qnorm(cdf_logitnormal(dados, par=p9ln$par))
# qqnorm(rln, main = "Logit-Normal")
# qqline(rln, col = "red", lwd = 2)

####  "Unit-Gompertz"
rug=qnorm(cdf_ugompertz(dados, par=p9ug$par))
# qqnorm(rug, main = "Unit-Gompertz")
# qqline(rug, col = "red", lwd = 2)

#### "Unit-Inverse Gaussian",
ruig=qnorm(cdf_uinvgauss(dados, par=p9uig$par))
# qqnorm(ruig, main = "Unit-Inverse Gaussian")
# qqline(ruig, col = "red", lwd = 2)

#### "Unit-Weibull"
ruw=qnorm(cdf_uweibull(dados, par=p9uw$par))
# qqnorm(ruw, main = "Unit-Weibull")
# qqline(ruw, col = "red", lwd = 2)


# UF
ruf=qnorm(cdf_uf(dados, par=p9uf$par))
# qqnorm(ruf, main = "Unit-Fréchet")
# qqline(ruf, col = "red", lwd = 2)


#### UBBS
rubbs=qnorm(cdf_ubbs(dados, par=p9ubbs$par))
# qqnorm(rubbs, main = "UBBS")
# qqline(rubbs, col = "red", lwd = 2)


#### Bimodal Beta
rbb=qnorm(cdf_bbeta(dados, par=p9bb$par))
# qqnorm(rbb, main = "Bimodal Beta")
# qqline(rbb, col = "red", lwd = 2)


# x11()
# par(mfrow=c(3,4))
# lista = list(rebb_finite, rtl, rul, rb, rkw, rln, rug, ruig, ruw, ruf, rubbs, rbb)
# for(i in 1:length(lista)){
#   qqenvelope(lista[[i]], ylab = "Randomized Quantile Residuals", main = "")
#   qqline(lista[[i]], col = 2)
# }


lista = list(rebb_finite, rtl, rul, rb, rkw, rln, rug, ruig, ruw, ruf, rubbs, rbb)
lista_nomes = c("EB", "Topp-Leone", "Unit-Lindley", "Beta", "Kumaraswamy", "Logit-Normal",
                "Unit-Gompertz", "Unit-Inverse Gaussian", "Unit-Weibull", "Unit-Fréchet",
                "UBBS", "Bimodal Beta")
x11()
par(mfrow=c(3,4))
for(i in 1:length(lista)){
  qqnorm(lista[[i]], main = lista_nomes[[i]],
         cex.main = 2,cex.lab = 1.5, cex.axis = 1.5)
  qqline(lista[[i]], col = "red", lwd = 2)
}


## replicar com boxplots
x11()
par(mfrow=c(3,4))
for(i in 1:length(lista)){
  boxplot(lista[[i]], main = lista_nomes[[i]], cex.main = 2)
}


###----------------------------------------------------------------
### Graficos complementares
###----------------------------------------------------------------
dados <- sort(d9)
n     <- length(dados)

# Parametros do modelo EB: c(alpha, beta, rho)
par_eb   <- c(p9b$par, p9ebb$rho_hat)
par_uig  <- p9uig$par
par_uf   <- p9uf$par
par_bb   <- p9bb$par


### --------------------------------------------------------------------------
### 1) CDFs ajustadas avaliadas nos dados
###    (o EB e reintegrado uma unica vez; os demais sao chamadas diretas)
### --------------------------------------------------------------------------

# EB: integracao numerica (gargalo, roda uma vez so)
#cdf_eb_dados  <- peb_integral(dados, par = par_eb)
cdf_eb_dados  <- F_Z(dados, par = par_eb)

cdf_uig_dados <- cdf_uinvgauss(par_uig, dados)
cdf_uf_dados  <- cdf_uf(par_uf,  dados)
cdf_bb_dados  <- cdf_bbeta(par_bb, dados)

# Lista organizada para reutilizacao
cdf_list <- list(
  "EB"                    = cdf_eb_dados,
  "Unit-Inverse Gaussian" = cdf_uig_dados,
  "Unit-Frechet"          = cdf_uf_dados,
  "Bimodal Beta"          = cdf_bb_dados
)

cores <- c("EB"                    = "red",
           "Unit-Inverse Gaussian" = "blue",
           "Unit-Frechet"          = "green",
           "Bimodal Beta"          = "purple")


### --------------------------------------------------------------------------
### 2) ECDF empirica
### --------------------------------------------------------------------------

Fn_func <- ecdf(dados)
Fn      <- Fn_func(dados)   # ECDF avaliada nos proprios dados


### ==========================================================================
### ANALISE 1 - GRAFICO DE DIFERENCA DE ECDF
### D(x) = F_empirica(x) - F_ajustada(x)
### Curvas proximas de zero => melhor ajuste. Mais informativo que sobrepor
### CDFs, pois amplifica os desvios que o olho nao capta na escala [0,1].
### ==========================================================================

# Matriz de diferencas (uma coluna por modelo)
dif_ecdf <- sapply(cdf_list, function(Fhat) Fn - Fhat)

# Estatistica-resumo: distancia integrada L1 e L2, e sup (tipo-KS),
# usadas como MEDIDAS DE DISTANCIA (nao como teste de hipotese).
resumo_dist <- data.frame(
  Modelo       = names(cdf_list),
  L1_integrada = apply(dif_ecdf, 2, function(d) mean(abs(d))),
  L2_integrada = apply(dif_ecdf, 2, function(d) sqrt(mean(d^2))),
  Sup_KS       = apply(dif_ecdf, 2, function(d) max(abs(d))),
  row.names    = NULL
)
cat("\n--- Distancias entre ECDF e CDF ajustada (menor = melhor) ---\n")
print(resumo_dist, digits = 4)

library(xtable)
xtable(resumo_dist, digits = 4, caption = "ECDF-difference")

# Grafico
x11(); par(mar = c(4.5, 4.5, 2, 1))
matplot(dados, dif_ecdf, type = "l", lty = 1, lwd = 2,
        col = cores[names(cdf_list)],
        xlab = "z", ylab = expression(F[n](z) - hat(F)(z)),
        main = "",
        cex.lab = 1.3, cex.main = 1.4)
abline(h = 0, col = "gray40", lty = 2)
legend("topright", legend = names(cdf_list),
       col = cores[names(cdf_list)], lty = 1, lwd = 2, cex = 1.1, bty = "n")


### ==========================================================================
### ANALISE 2 - PIT (PROBABILITY INTEGRAL TRANSFORM)
### Se o modelo esta correto, u_i = F(x_i; theta_hat) ~ Uniforme(0,1).
### Desvios da uniformidade => ma especificacao. Histograma + linha de
### referencia da densidade uniforme (= 1).
### ==========================================================================

x11(); par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (nome in names(cdf_list)) {
  u <- cdf_list[[nome]]
  u <- u[is.finite(u) & u > 0 & u < 1]   # descarta saturacoes eventuais
  hist(u, breaks = 20, probability = TRUE,
       col = "gray85", border = "white",
       xlim = c(0, 1), ylim = c(0, 2),
       main = nome, xlab = "PIT", cex.main = 1.4)
  abline(h = 1, col = "red", lwd = 2, lty = 1)  # uniforme ideal
}

# Resumo quantitativo do PIT: desvio da uniformidade via distancia
# de Kolmogorov contra a U(0,1) (reportada como MEDIDA, nao p-valor).
resumo_pit <- data.frame(
  Modelo   = names(cdf_list),
  KS_vs_Unif = sapply(cdf_list, function(u) {
    u <- sort(u[is.finite(u) & u > 0 & u < 1])
    m <- length(u)
    max(abs(u - (seq_len(m) - 0.5)/m))   # dist. sup a U(0,1)
  }),
  row.names = NULL
)
print(resumo_pit, digits = 4)

xtable(resumo_pit, digits = 4)


tabela_diagnosticos <- merge(resumo_dist, resumo_pit, by = "Modelo")
print(tabela_diagnosticos, digits = 4, row.names = FALSE)

xtable(tabela_diagnosticos, digits = 2, caption = "ECDF-difference")

# Salvar em disco (ajuste o caminho conforme necessario)
# write.csv2(tabela_diagnosticos, "tabela_diagnosticos_revisor.csv", row.names = FALSE)

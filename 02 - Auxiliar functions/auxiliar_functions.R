library(gsl)
library(MASS)

###-----------------------------------------------
### Funcoes auxiliares
###-----------------------------------------------

## Beta function
B=function(a,b){gamma(a)*gamma(b)/gamma(a+b)}

#incomplete beta function
IB <- function(x, a, b) {
  return(pbeta(x, a, b) * beta(a, b))
}
## Regularized incomplete beta function
# library(zipfR)
# Ibeta(x, a, b)
rIB <- function(x, a, b) {
  return(pbeta(x, a, b))  
}

# Ibeta(0.5, 2, 3)
# IB(0.5, 2, 3)  # Returns I_x(2,3)
# rIB(0.5, 2, 3)  # Returns I_x(2,3)

# Example usage:
#rIB(0.5, 2, 3)  # Returns I_x(2,3)

##2F1 function (hypergeometric function)
#library(gsl)
# hyperg_2F1(a, b, c, x)
# hyperg_2F1(1, 1, 1, 1);hyperg_2F1(1, 1, 1, 0);hyperg_2F1(1, 1, 1, 0.5)



# The Appell function F2

# Appell F2 function
AppellF2 <- function(a, b1, b2, c1, c2, x, y, tol = 1e-12, max_terms = 200) {
  poch <- function(q, n) gamma(q + n) / gamma(q)
  
  total <- 0
  for (m in 0:max_terms) {
    for (n in 0:max_terms) {
      term <- poch(a, m + n) * poch(b1, m) * poch(b2, n) /
        (poch(c1, m) * poch(c2, n) * factorial(m) * factorial(n)) *
        (x^m) * (y^n)
      
      # Check for numerical problems
      if (!is.finite(term)) next
      
      total <- total + term
      
      # Convergence check only if finite
      if (abs(term) < tol) break
    }
  }
  return(total)
}


#
# Test example
# AppellF2(1, 1, 1, 2, 2, 0.2, 0.3)
# 
# # Example: F2(1, 1, 1, 2, 2; 0.2, 0.3)
# F2(1, 1, 1, 2, 2, 0.2, 0.3)
# ## Resultado do Wolfram
# #AppellF2[2, 1, 3, 4, 5, 0.6, 0.3]
# #2.61933
# #F2_serie(0.6, 0.3, 2,1,3,4,5)#2.619318
# AppellF2(2,1,3,4,5, 0.6, 0.3)


## Condicoes de convergencia dessa serie: |x|+|y|<1
## Caso outra condicao seja requerida, representacoes integrais serao necessárias.
## No caso da extended beta distribution, essa condicao eh atendida para CDF e PDF




# -----------------------------------------------------------------------------
# CDF
# -----------------------------------------------------------------------------
F_Z <- function(par, x) {
  alpha <- par[1]; beta <- par[2]; rho <- par[3]
  
  integrand <- function(u, s) {
    G_a <- pgamma(u,     shape = alpha)
    G_b <- pgamma(s * u, shape = beta)
    u^(alpha - 1) * exp(-u) *
      (G_b + rho * (2 * G_a - 1) * G_b * (G_b - 1))
  }
  
  sapply(x, function(zi) {
    if (zi <= 0) return(0)
    if (zi >= 1) return(1)
    s   <- 1 / zi - 1
    val <- tryCatch(
      integrate(integrand, 0, Inf, s = s,
                rel.tol = 1e-8, abs.tol = 1e-12,
                subdivisions = 1000L)$value,
      error = function(e) NA_real_
    )
    if (is.na(val)) return(NA_real_)
    1 - val / gamma(alpha)
  })
}

# -----------------------------------------------------------------------------
# Densidade
# f(z) = (1+rho)*K1 - rho*K2 - rho*K3 + rho*K4
# -----------------------------------------------------------------------------
f_Z <- function(par, x) {
  alpha <- par[1]; beta <- par[2]; rho <- par[3]
  
  Ki <- function(s, tipo) {
    C   <- (1 + s)^2 * s^(beta - 1) / (gamma(alpha) * gamma(beta))
    nuc <- function(u) u^(alpha + beta - 1) * exp(-(1 + s) * u)
    
    integrando <- switch(tipo,
                         K1 = function(u) nuc(u),
                         K2 = function(u) nuc(u) * pgamma(u,     shape = alpha),
                         K3 = function(u) nuc(u) * pgamma(s * u, shape = beta),
                         K4 = function(u) nuc(u) * pgamma(u, shape = alpha) *
                           pgamma(s * u, shape = beta)
    )
    fat <- c(K1 = 1, K2 = 2, K3 = 2, K4 = 4)[tipo]
    val <- tryCatch(
      integrate(integrando, 0, Inf,
                rel.tol = 1e-8, abs.tol = 1e-12,
                subdivisions = 1000L)$value,
      error = function(e) NA_real_
    )
    C * fat * val
  }
  
  sapply(x, function(zi) {
    if (zi <= 0 || zi >= 1) return(0)
    s  <- 1 / zi - 1
    k1 <- Ki(s, "K1"); k2 <- Ki(s, "K2")
    k3 <- Ki(s, "K3"); k4 <- Ki(s, "K4")
    val <- (1 + rho) * k1 - rho * k2 - rho * k3 + rho * k4
    if (!is.finite(val)) NA_real_ else val
  })
}

# =============================================================================
# Verificações
# =============================================================================
# par0   <- c(alpha = 2, beta = 3, rho = 0.5)
# z_test <- c(0.1, 0.3, 0.5, 0.7, 0.9)
# 
# cat("=== CDF ===\n")
# print(round(F_Z(par0, z_test), 6))
# 
# cat("\n=== Densidade ===\n")
# print(round(f_Z(par0, z_test), 6))
# 
# cat("\n=== Integral da densidade (deve ser 1) ===\n")
# for (r in c(-0.9, -0.5, 0, 0.5, 0.9)) {
#   val <- tryCatch(
#     integrate(f_Z, 0.01, 0.99, par = c(2, 3, r),
#               rel.tol = 1e-6, subdivisions = 500L)$value,
#     error = function(e) NA_real_
#   )
#   cat(sprintf("rho = %+.1f  =>  %.8f\n", r, val))
# }




deb_fast <- function(par, x) {
  a <- par[1]  # alpha
  b <- par[2]  # beta
  r <- par[3]  # rho
  
  # Pré-calcular valores constantes
  ga <- gamma(a)
  gb <- gamma(b)
  gab <- B(a, b)
  ga2b <- gamma(2 * a + b)
  gab2 <- gamma(a + 2 * b)
  ga2b2 <- gamma(2 * a + 2 * b)
  
  # Pré-calcular coeficientes fixos
  n21 <- 2 * r * ga2b / (a * ga^2 * gb)
  n31 <- 2 * r * gab2 / (b * ga * gb^2)
  n41 <- (4^(1 - a - b)) * r * ga2b2 / (a * b * ga^2 * gb^2)
  
  # Filtrar x válidos
  mask <- x > 0 & x < 1
  out <- numeric(length(x))
  if (!any(mask)) return(out)
  
  z <- x[mask]
  
  # Pré-computar partes fáceis
  z1 <- z^(a - 1)
  z2 <- (1 - z)^(b - 1)
  B_ab <- gab
  
  n1 <- (1 + r) * z1 * z2 / B_ab
  
  # Parte com 2F1
  n22 <- hyperg_2F1(1, 2 * a + b, 1 + a, z / (1 + z))
  n23 <- z^(2 * a - 1) * (1 - z)^(b - 1) / (1 + z)^(2 * a + b)
  n2 <- n21 * n22 * n23
  
  n32 <- hyperg_2F1(1, a + 2 * b, 1 + b, (1 - z) / (2 - z))
  n33 <- z^(a - 1) * (1 - z)^(2 * b - 1) / (2 - z)^(a + 2 * b)
  n3 <- n31 * n32 * n33
  
  # Appell F2 — pode ser MUITO lento se não otimizado
  n42 <- vapply(z, function(zz) {
    AppellF2(2 * a + 2 * b, 1, 1, a + 1, b + 1, zz / 2, (1 - zz) / 2)
  }, numeric(1))
  n43 <- z^(2 * a - 1) * (1 - z)^(2 * b - 1)
  n4 <- n41 * n42 * n43
  
  out[mask] <- n1 - n2 - n3 + n4
  return(out)
}

#deb_fast(seq(0,1,0.1), c(2.5, 2.5, -0.5))
#deb_fast_r(seq(0,1,0.1), c(2.5, 2.5, -0.5))

###-----------------------------------------------

peb_integral <- function(x, par) {
  # Garante vetor de saída
  Fvals <- numeric(length(x))
  
  # Garante ordem crescente (para integração cumulativa)
  ord <- order(x)
  x_ord <- x[ord]
  
  # Loop em cada ponto (pode ser vetorizado depois)
  for (i in seq_along(x_ord)) {
    z <- x_ord[i]
    
    if (z <= 0) {
      Fvals[i] <- 0
    } else if (z >= 1) {
      Fvals[i] <- 1
    } else {
      # Integração numérica de deb() de 0 até z
      val <- tryCatch(
        integrate(function(t) deb_fast(t, par),
                  lower = 0, upper = z,
                  rel.tol = 1e-6, abs.tol = 1e-8)$value,
        error = function(e) NA
      )
      Fvals[i] <- val
    }
  }
  
  # Normalização: F(1) deve ser 1
  if (max(Fvals, na.rm = TRUE) > 0)
    Fvals <- Fvals / max(Fvals, na.rm = TRUE)
  
  # Reordena de volta
  Fout <- numeric(length(x))
  Fout[ord] <- Fvals
  
  # Força limites [0,1] e monotonicidade
  Fout <- pmin(pmax(Fout, 0), 1)
  Fout <- cummax(Fout)
  
  return(Fout)
}

# peb_integral(seq(0,1,0.1),c(2.5, 2.5, -0.5))



######################################
### random sample generator
######################################


# Metodo de otimizacao (Gradient descent)
library(numDeriv)

library(MASS)
#?fitdistr
###-----------------------------------------------
## Gerador de amostras (baseado em (X,Y))
###-----------------------------------------------

# Função para calcular a inversa da CDF da distribuição gamma
gamma_inv_cdf <- function(p, shape, scale) {
  qgamma(p, shape = shape, scale = scale)
}



reb2 <- function(n, par){
  alpha=par[1]#alpha
  beta=par[2]#beta
  rho=par[3]
  
  
  
  # Vetores para armazenar amostras
  X_samples <- numeric(n)
  Y_samples <- numeric(n)
  Z_samples <- numeric(n)
  
  for (i in 1:n) {
    # Passo 1: Simular U2 ~ U(0,1)
    U2 <- runif(1)
    
    # Passo 2: Definir a CDF condicional FU1|U2
    FU1_U2 <- function(u1) {
      u1 * (1 + rho * (1 - u1) * (1 - 2 * U2))
      #u1 * (1 + 2*rho * (1 - u1) * (1 - 2 * U2))# teste pois no anterior todos os estimadores estao obtendo rho/2
    }
    
    # # Passo 3: Calcular inversa generalizada numericamente
    inverse_FU1_U2 <- function(v) {
      uniroot(function(u1) FU1_U2(u1) - v, interval = c(0, 1))$root
    }
    
    # Passo 4: Gerar V ~ U(0,1)
    V <- runif(1)
    
    # Passo 5: Definir U1
    U1 <- inverse_FU1_U2(V)
    
    # Passo 6: Obter (X, Y) das distribuições gama
    X <- gamma_inv_cdf(U1, shape = alpha, scale = 1)  # Escala assumida como 1
    Y <- gamma_inv_cdf(U2, shape = beta, scale = 1)   # Escala assumida como 1
    
    # Passo 7: Calcular Z
    Z_samples[i] <- X / (X + Y)
    X_samples[i] <- X
    Y_samples[i] <- Y
  }
  
  return(list(Z=Z_samples, X=X_samples, Y=Y_samples))
}

# # Exemplo de uso
# set.seed(123)
# reb2(1000, c(2, 3, 0.5))
# x11()
# # Visualizar a amostra
# hist(reb2(1000, c(2, 3, 0.5))$Z, breaks = 30, main = "Distribuição de Z", col = "lightblue", probability = TRUE)
# hist(reb2(1000, c(1.5, 2.5, -0.9))$Z)

###-----------------------------------------------
## Estimação de parâmetros
###-----------------------------------------------

### caso de dados bivariados correlacionados (X,Y) --> Z=X/(X+Y)


#library(MASS)
##funcao de densidade conjunta
f_XY = function(x,y, r, par_est){
  alpha_est = par_est[1]
  beta_est = par_est[2]
  
  out = NULL
  for(i in 1:length(x)){
    fx = dgamma(x[i], shape = alpha_est, rate = 1)
    fy = dgamma(y[i], shape = beta_est, rate = 1)
    Fx=pgamma(x[i], shape = alpha_est, rate = 1)
    Fy=pgamma(y[i], shape = beta_est, rate = 1)
    out = c(out, fx*fy*(1+r*(2*Fx-1)*(2*Fy-1)))
  }
  return(out)
  
}


### corrigindo o espaco de rho [-1,1] via projection step
estimador_r_bivariado <- function(data, w0 = NULL,
                                  max_iter = 100, step_size = 0.01) {
  x <- data$X
  y <- data$Y
  
  # estimar parâmetros marginais
  alpha_est <- as.numeric(fitdistr(x, "gamma")$estimate[1])
  beta_est  <- as.numeric(fitdistr(y, "gamma")$estimate[1])
  
  # inicialização
  #w <- w0
  w <- ifelse(is.null(w0), 0, w0)
  ws <- numeric(max_iter)
  motivo_parada <- "convergência normal"
  
  for (i in seq_len(max_iter)) {
    # tentativa segura de calcular gradiente
    v <- tryCatch(
      grad(function(r) -sum(log(f_XY(x, y, r, c(alpha_est, beta_est)))), w),
      error = function(e) NA_real_
    )
    
    # se o gradiente é inválido, interrompe
    if (!is.finite(v)) {
      motivo_parada <- paste("gradiente não finito na iteração", i)
      break
    }
    
    # atualização
    w_new <- w - step_size * v
    
    # projeção no intervalo [-1,1]
    if (w_new < -1) {
      w_new <- -1
    } else if (w_new > 1) {
      w_new <- 1
    }
    
    # verificação de validade
    if (!is.finite(w_new)) {
      motivo_parada <- paste("w não finito na iteração", i)
      break
    }
    
    # salvar histórico
    ws[i] <- w_new
    w <- w_new
  }
  
  iter_realizadas <- i
  
  # resultado final
  return(list(
    alpha = alpha_est,
    beta = beta_est,
    w_final = w,
    iteracoes = iter_realizadas,
    motivo_parada = motivo_parada,
    w_hist = ws[1:iter_realizadas]
  ))
}


# set.seed(1)
# estimador_r_bivariado(reb2(100, c(1.5, 2.5, -0.7)))
# estimador_r_bivariado(reb2(1000, c(1.5, 2.5, -0.7)))

###-----------------------------------------------
### Estimacao em caso de dados unitarios
###-----------------------------------------------


#### Estimador simples para rho (apos estimar alphas e mesma escala ajustada)
estimador_rho_gamma3 <- function(x, y, alpha1, alpha2) {
  stopifnot(length(x) == length(y))
  
  # CDFs teóricas
  Gx <- pgamma(x, shape = alpha1, scale = 1)
  Hy <- pgamma(y, shape = alpha2, scale = 1)
  
  # Corr de Spearman-like via CDFs
  rho_hat <- cor(Gx, Hy)
  
  # Beta_x e Beta_y (teóricos)
  G <- function(x) pgamma(x, shape = alpha1, scale = 1)
  H <- function(y) pgamma(y, shape = alpha2, scale = 1)
  integrand_G <- function(x) G(x) * (1 - G(x))
  integrand_H <- function(y) H(y) * (1 - H(y))
  beta_x <- integrate(integrand_G, lower = 0, upper = Inf)$value
  beta_y <- integrate(integrand_H, lower = 0, upper = Inf)$value
  
  # Estatísticas só para referência
  sx <- sd(x)
  sy <- sd(y)
  
  return(list(
    rho_hat = rho_hat,
    cor_xy = cor(x, y),
    beta_x = beta_x,
    beta_y = beta_y,
    sx = sx,
    sy = sy
  ))
}



ebb_fit <- function(data, par1){
  loglik_r <- function(r){
    vals <- deb_fast(data, par=c(par1, r))
    vals <- pmax(vals, 1e-300)
    sum(log(vals))
  }
  
  res <- optimize(loglik_r, interval = c(-0.999, 0.999), maximum = TRUE)
  list(r_opt = res$maximum, llf = res$objective)
  return(list(par=c(par1,res$maximum), llf_max=res$objective, 
              AIC=aic(res$objective, 3),
              BIC=bic(res$objective, 3, length(data))))
}


###-----------------------------------------------
## Gerador de amostras (baseado em (X,Y))
###-----------------------------------------------

# Função para calcular a inversa da CDF da distribuição gamma
gamma_inv_cdf <- function(p, shape, scale) {
  qgamma(p, shape = shape, scale = scale)
}

### -----------------------------------------
### Ajuste dos modelos concorrentes

# library(gsl)
# library(MASS)
# library(latex2exp)
# library(AdequacyModel)
# library(goftest)
# 
# #Carrega os pacotes necessários para realizar o paralelismo
# library(foreach)
# library(doParallel)
### ----------------------------------------------
# Ajuste dos modelos concorrentes
### ----------------------------------------------
# 
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


### ---------------------------------------------
### Ajustes dos modelos beta e kumaraswamy
### ---------------------------------------------
library(AdequacyModel)
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




#################################################
# 1) BETA
#################################################

pdf_beta <- function(par, x){
  a = par[1]
  b = par[2]
  
  fx = dbeta(x, a, b)
  
  return(
    ifelse(x > 0 & x < 1,
           fx,
           0)
  )
}

cdf_beta <- function(par, x){
  a = par[1]
  b = par[2]
  
  Fx = pbeta(x, a, b)
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}

#library(AdequacyModel)
fit_beta=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_beta, cdf = cdf_beta,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w)),
              result_1))
}

#################################################
# 2) KUMARASWAMY
#################################################

pdf_kw <- function(par, x){
  a = par[1]
  b = par[2]
  
  fx = a*b*(x^(a-1))*((1-x^a)^(b-1))
  
  return(
    ifelse(x > 0 & x < 1,
           fx,
           0)
  )
}

cdf_kw <- function(par, x){
  a = par[1]
  b = par[2]
  
  Fx = 1 - (1 - x^a)^b
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}


fit_kw=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_kw, cdf = cdf_kw,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}



#################################################
# 3) LOGIT-NORMAL
#################################################

pdf_logitnormal <- function(par, x){
  mu    = par[1]
  sigma = par[2]
  
  num = exp(-((log(x/(1-x)) - mu)^2)/(2*sigma^2))
  den = sigma*sqrt(2*pi)*x*(1-x)
  
  return(
    ifelse(x > 0 & x < 1,
           num/den,
           0)
  )
}

cdf_logitnormal <- function(par, x){
  mu    = par[1]
  sigma = par[2]
  
  z = (log(x/(1-x)) - mu)/sigma
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  pnorm(z)))
  )
}


fit_logitnormal=function(w, start_values=c(1,1), lim_inf= c(-30,0),
                         lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_logitnormal, cdf = cdf_logitnormal,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}

#################################################
# 4) TRUNCATED NORMAL
#################################################

pdf_tnormal <- function(par, x){
  mu    = par[1]
  sigma = par[2]
  
  A = pnorm((1-mu)/sigma) - pnorm((-mu)/sigma)
  
  num = dnorm((x-mu)/sigma)
  den = sigma*A
  
  return(
    ifelse(x > 0 & x < 1,
           num/den,
           0)
  )
}

cdf_tnormal <- function(par, x){
  mu    = par[1]
  sigma = par[2]
  
  A = pnorm((1-mu)/sigma) - pnorm((-mu)/sigma)
  
  Fx = (pnorm((x-mu)/sigma) - pnorm((-mu)/sigma))/A
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}


fit_tnormal=function(w, start_values=c(1,1), lim_inf= c(-30,0), lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_tnormal, cdf = cdf_tnormal,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}


#################################################
# 5) TOPP-LEONE
#################################################

pdf_toppleone <- function(par, x){
  a = par[1]
  
  fx = 2*a*(x^(a-1))*((2-x)^a)*(1-x)
  
  return(
    ifelse(x > 0 & x < 1,
           fx,
           0)
  )
}

cdf_toppleone <- function(par, x){
  a = par[1]
  
  Fx = (x*(2-x))^a
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}

fit_toppleone=function(w, start_values=c(1), lim_inf= c(0), lim_sup=c(30)){
  result_1 = goodness.fit(pdf = pdf_toppleone, cdf = cdf_toppleone,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 1),
              BIC=bic(-result_1$Value, 1, length(w))))
}


#################################################
# 6) UNIT-LINDLEY
#################################################

pdf_ulindley <- function(par, x){
  theta = par[1]
  
  fx = (theta^2/(1+theta)) *
    (1-x)^(-3) *
    exp(-(theta*x)/(1-x))
  
  return(
    ifelse(x > 0 & x < 1,
           fx,
           0)
  )
}

cdf_ulindley <- function(par, x){
  theta = par[1]
  
  Fx = 1 - (1 + (theta*x)/(1-x)) *
    exp(-(theta*x)/(1-x))
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}

fit_ulindley=function(w, start_values=c(1), lim_inf= c(0), lim_sup=c(30)){
  result_1 = goodness.fit(pdf = pdf_ulindley, cdf = cdf_ulindley,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1), mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 1),
              BIC=bic(-result_1$Value, 1, length(w))))
}

#################################################
# 7) UNIT-GOMPERTZ
#################################################

# pdf_ugompertz <- function(par, x){
#   alpha = par[1]
#   beta  = par[2]
#   
#   fx = alpha*beta*
#     x^(-(beta+1)) *
#     exp(-alpha*(x^(-beta)-1))
#   
#   return(
#     ifelse(x > 0 & x < 1,
#            fx,
#            0)
#   )
# }
# 
# cdf_ugompertz <- function(par, x){
#   alpha = par[1]
#   beta  = par[2]
#   
#   Fx = exp(-alpha*((x^(-beta))-1))
#   
#   return(
#     ifelse(x <= 0, 0,
#            ifelse(x >= 1, 1,
#                   Fx))
#   )
# }

pdf_ugompertz <- function(par, x){
  a = par[1];b=par[2]
  return(ifelse(x>0 & x<1,  a*b*exp(-a*((x^(-b)) - 1))/(x^(1+b)), 0))
}

cdf_ugompertz <- function(par, x){
  a = par[1];b=par[2]
  return(ifelse(x<=0, 0, ifelse(x>=1, 1, exp(-a*((x^(-b)) - 1)))))
}



fit_ugompertz=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_ugompertz, cdf = cdf_ugompertz,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}

#################################################
# 8) UNIT-INVERSE GAUSSIAN
#################################################

pdf_uinvgauss <- function(par, x){
  mu     = par[1]
  lambda = par[2]
  
  c1 = sqrt(lambda/(2*pi))
  
  fx = c1 *
    (x*(1-x))^(-3/2) *
    exp(-(lambda*(x - mu*(1-x))^2)/
          (2*(mu^2)*x*(1-x)))
  
  return(
    ifelse(x > 0 & x < 1,
           fx,
           0)
  )
}

cdf_uinvgauss <- function(par, x){
  mu     = par[1]
  lambda = par[2]
  
  z1 = sqrt((lambda*(1-x))/x) *
    (x/(mu*(1-x)) - 1)
  
  z2 = -sqrt((lambda*(1-x))/x) *
    (x/(mu*(1-x)) + 1)
  
  Fx = pnorm(z1) +
    exp(2*lambda/mu) * pnorm(z2)
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}

fit_uinvgauss=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_uinvgauss, cdf = cdf_uinvgauss,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}


#################################################
# 8) Bimodal Beta
#################################################
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
  return(ifelse(x>0 & x<1, c1*(x^(a-1))*(1-x)^(b-1), 0))
}
# cdf_bbeta <- function(par,x){
#   a = par[1]
#   b = par[2]
#   r=par[3]
#   d=par[4]
#   
#   z=1+r-2*d*(a/(a+b))+(d^2)*(a*(a+1))/((a+b)*(a+b+1))
#   
#   c1=(1+r)*Rbeta(x,a,b) - (2*d*Ibeta(x,a,b)/beta(a,b)) + (d^2)*Ibeta(x,a+2,b)/beta(a,b)
#   
#   return(  ifelse(x<0, 0,
#                   ifelse(x>1, 1,c1/z)))
# }
cdf_bbeta <- function(par, x){
  
  a = par[1]
  b = par[2]
  r = par[3]
  d = par[4]
  
  z = 1 + r -
    2*d*(a/(a+b)) +
    d^2*(a*(a+1))/((a+b)*(a+b+1))
  
  out = numeric(length(x))
  
  out[x <= 0] = 0
  out[x >= 1] = 1
  
  ind = (x > 0 & x < 1)
  
  if(any(ind)){
    
    xx = x[ind]
    
    c1 =
      (1+r)*Rbeta(xx,a,b) -
      (2*d*Ibeta(xx,a+1,b)/beta(a,b)) +
      (d^2*Ibeta(xx,a+2,b)/beta(a,b))
    
    out[ind] = c1/z
  }
  
  return(out)
}
#cdf_bbeta(c(1,1,0.5,1), runif(100))

fit_bbeta=function(w, start_values=c(1,1,1,0),
                   lim_inf= c(0,0,0,-30), lim_sup=c(30,30,30,30)){
  result_1 = goodness.fit(pdf = pdf_bbeta, cdf = cdf_bbeta,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 4),
              BIC=bic(-result_1$Value, 4, length(w))))
}


#################################################
# 9) Unit-Weibull
#################################################

pdf_uweibull <- function(par, x){
  a = par[1]
  b = par[2]
  
  fx = (1/x)*a*b*((-log(x))^(b-1))*exp(-a*((-log(x))^b))
  
  return(
    ifelse(x > 0 & x < 1,
           fx,
           0)
  )
}

cdf_uweibull <- function(par, x){
  a = par[1]
  b = par[2]
  
  Fx = exp(-a*((-log(x))^b))
  
  return(
    ifelse(x <= 0, 0,
           ifelse(x >= 1, 1,
                  Fx))
  )
}

fit_uweibull=function(w, start_values=c(1,1), lim_inf= c(0,0), lim_sup=c(30,30)){
  result_1 = goodness.fit(pdf = pdf_uweibull, cdf = cdf_uweibull,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 2),
              BIC=bic(-result_1$Value, 2, length(w))))
}

#################################################
# 9) UBBS
#################################################

pdf_ubbs <- function(par,x){
  a = par[1]
  b = par[2]
  d = par[3]
  
  s=(1/a)*(sqrt(-log(x)/b)-sqrt(-b/log(x)))
  
  n1=((-b/log(x))^(1/2))+((-b/log(x))^(3/2))
  q1=4*x*a*b*pnorm(-d)
  n2=dnorm(abs(s)+d)
  
  out=n1*n2/q1
  
  out2=ifelse(x>=1, 1, ifelse(x<=0, 0, out))
  return(out2)
}
#pdf_ubbs(c(1,1,1), runif(10))## vetorizada!!

cdf_ubbs <- function(par,x){
  a = par[1]
  b = par[2]
  d = par[3]
  
  out = rep(0, length(x))
  
  ind=which(x>0 & x<1)
  
  s=(1/a)*(sqrt(-log(x[ind])/b)-sqrt(-b/log(x[ind])))
  out=ifelse(x[ind]<=exp(-b),
             (1-pnorm(s+d))/(2*(1-pnorm(d))),
             1-((pnorm(s-d))/(2*(1-pnorm(d)))))
  
  out2=ifelse(x>=1, 1, ifelse(x<=0, 0, out))
  
  return(out2)
  
}
#cdf_ubbs(c(1,1,1), runif(10))## vetorizada!!

fit_ubbs=function(w, start_values=c(1,1,0.1), lim_inf= c(0,0,-3), lim_sup=c(30,30,3)){
  result_1 = goodness.fit(pdf = pdf_ubbs, cdf = cdf_ubbs,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 3),
              BIC=bic(-result_1$Value, 4, length(w))))
}


#################################################
# 10) unit Frechet distribution
#################################################
#par: (sigma, alpha, rho), s.t. sigma>0, alpha>0 and 0<=rho<=1
#w: data (values in (0,1))
pdf_uf<- function(par,w){
  s1=par[1]#sigma
  a=par[2]#alpha
  rho=par[3]#rho
  
  s=(w/(1-w))
  
  c1=a/(s1^a)* (s^{a-1})* ((s+1)^2)
  c2=2*((((s^a)/(s1^a))+1)^2) - rho*(((s^a)/(s1^a))^2 +1)
  d1=((((s^a)/(s1^a)+1)^2) -rho*(s^a)/(s1^a)  )^2
  d2=(((s^a)/(s1^a))+1)^2
  
  out=c1*((c2/d1)-1/d2)
  return(out)
}


#par: (sigma, alpha, rho), s.t. sigma>0, alpha>0 and 0<=rho<=1
#w: data (values in (0,1))
cdf_uf<- function(par, x){
  w=x
  
  s1=par[1]#sigma
  a=par[2]#alpha
  rho=par[3]#rho
  
  s=(w/(1-w))
  
  c1=(s^a)/(s1^a)
  c2=((((s^a)/(s1^a))+1)^2)-rho 
  d1=((((s^a)/(s1^a)+1)^2) -rho*(s^a)/(s1^a)  )
  d2=(((s^a)/(s1^a))+1)
  
  out=c1*c2/(d1*d2)
  
  out=ifelse(w<=0, 0, ifelse(w>=1,1, out))
  return(out)
  
}

fit_uf=function(w, start_values=c(1,1), lim_inf= c(0,0,0), lim_sup=c(30,30,1)){
  result_1 = goodness.fit(pdf = pdf_uf, cdf = cdf_uf,
                          starts = start_values, data = w, method = "PSO",
                          domain = c(0,1),mle = NULL, lim_inf =lim_inf,
                          lim_sup = lim_sup, S = 250, prop=0.1, N=50)
  return(list(par=result_1$mle, llf_max=-result_1$Value, 
              AIC=aic(-result_1$Value, 3),
              BIC=bic(-result_1$Value, 3, length(w))))
}



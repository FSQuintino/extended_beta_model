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


#PDF
deb_fast <- function(x, par) {
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


# CDF via numerical integral
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



######################################
### random sample generator
######################################

###-----------------------------------------------
## random sample generator (based on (X,Y))
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

# # Example
# set.seed(123)
# reb2(1000, c(2, 3, 0.5))
# x11()
# # plots
# hist(reb2(1000, c(2, 3, 0.5))$Z, breaks = 30, col = "lightblue", probability = TRUE)
# hist(reb2(1000, c(1.5, 2.5, -0.9))$Z)

###-----------------------------------------------
## Parameter estimation
###-----------------------------------------------

### bivariate case: (X,Y) --> Z=X/(X+Y)


#library(MASS)
##joint pdf
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



### estimator
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


### MLE for DEB model
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

### -------------------------------------------
## Metodos concorrentes para ajuste do Modelo Beta
### -------------------------------------------

#estimador via minimizacao da diferenca entre densidade kernel e PDF Beta
estim_beta_kernel <- function(x){
  
  # densidade kernel em uma grade de 200 pontos
  dens_k <- density(x, n = 200)
  t_grid <- dens_k$x
  f_grid <- dens_k$y
  
  # função de perda a minimizar
  loss_fun <- function(par){
    a <- par[1]
    b <- par[2]
    
    # penaliza parâmetros inválidos
    if(a <= 0 || b <= 0) return(1e10)
    
    f_beta <- dbeta(t_grid, a, b)
    sum( (f_grid - f_beta)^2 )
  }
  
  # chute inicial
  init <- c(1, 1)
  
  # otimização
  opt <- optim(
    par = init,
    fn = loss_fun,
    method = "L-BFGS-B",
    lower = c(0.01, 0.01),
    upper = c(50, 50)
  )
  
  return(list(alpha = opt$par[1],
              beta  = opt$par[2],
              value = opt$value))
}


# set.seed(1)
# estim_beta_kernel(rbeta(1000, 3, 10))
# estim_beta_kernel(rbeta(30, 3, 10))

# estimador via minimizacao da diferenca entre ECDF e CDF Beta
estim_beta_ecdf <- function(x){
  
  # ECDF
  Fn <- ecdf(x)
  
  # grade de pontos (pode ajustar)
  t_grid <- seq(0.001, 0.999, length.out = 300)
  F_emp  <- Fn(t_grid)
  
  # função de perda
  loss_fun <- function(par){
    a <- par[1]
    b <- par[2]
    
    # penaliza parâmetros inválidos
    if(a <= 0 || b <= 0) return(1e12)
    
    F_beta <- pbeta(t_grid, a, b)
    sum( (F_emp - F_beta)^2 )
  }
  
  # chute inicial
  init <- c(1, 1)
  
  # otimização
  opt <- optim(
    par   = init,
    fn    = loss_fun,
    method = "L-BFGS-B",
    lower = c(0.01, 0.01),
    upper = c(50, 50)
  )
  
  return(list(alpha = opt$par[1],
              beta  = opt$par[2],
              value = opt$value))
}

# set.seed(1)
# estim_beta_ecdf(rbeta(30, 3, 10))
# estim_beta_ecdf(rbeta(1000, 3, 10))

# estimador Bayesiano MAP para parametros da Beta
beta_bayes_map <- function(x, a0 = 1, b0 = 1, c0 = 1, d0 = 1){
  
  logpost <- function(par){
    a <- par[1]
    b <- par[2]
    
    if(a <= 0 || b <= 0) return(-Inf)
    
    n <- length(x)
    
    # log-verossimilhança Beta
    ll <- sum( dbeta(x, a, b, log = TRUE) )
    
    # priors independentes: Gamma(a0,b0), Gamma(c0,d0)
    lp <- dgamma(a, shape = a0, rate = b0, log = TRUE) +
      dgamma(b, shape = c0, rate = d0, log = TRUE)
    
    return(-(ll + lp))  # negativo: optim minimiza
  }
  
  opt <- optim(
    par = c(1,1),
    fn  = logpost,
    method = "L-BFGS-B",
    lower = c(1e-4, 1e-4),
    upper = c(200, 200)
  )
  
  list(alpha = opt$par[1],
       beta  = opt$par[2],
       logpost = -opt$value)
}

# set.seed(1)
#  beta_bayes_map(rbeta(30, 3, 10))
#  beta_bayes_map(rbeta(1000, 3, 10))



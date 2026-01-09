###-----------------------------------------------
### auxiliar functions
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
library(gsl)
# hyperg_2F1(a, b, c, x)
hyperg_2F1(1, 1, 1, 1);hyperg_2F1(1, 1, 1, 0);hyperg_2F1(1, 1, 1, 0.5)



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


## Condicoes de convergencia dessa serie: |x|+|y|<1
## Caso outra condicao seja requerida, representacoes integrais serao necessárias.
## No caso da extended beta distribution, essa condicao eh atendida para CDF e PDF


###-----------------------------------------------
#PDF of extended beta distribution

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
#deb_fast(seq(0,1,0.1), c(2.5, 2.5, -0.5))






####----------------------------------------
#### Graficos das PDFs
####----------------------------------------
library(latex2exp)

#################################
### Vatiando alpha
x11()
par(mfrow=c(1,3)) 

x_points = seq(0,1, by=0.01)
#x_points=x_points[-c(1,length(x_points))]

par1=c(0.5, 1.5,0.5)
par2=c(1, 1.5,0.5)
par3=c(1.5, 1.5,0.5)
par4=c(2, 1.5,0.5)
par5=c(2.5, 1.5,0.5)


(g1=deb_fast(x_points, par1))
(g2=deb_fast(x_points, par2))
(g3=deb_fast(x_points, par3))
(g4=deb_fast(x_points, par4))
(g5=deb_fast(x_points, par5))


plot(x_points, g1,
     ylim = c(0,8),
     main=TeX("$beta=1.5, rho=0.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g2, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$alpha$"), 
       c(0.5,1,1.5,2, 2.5), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)


#################################
### Vatiando beta
# x11()
# par(mfrow=c(1,2)) 

# x_points = seq(0,1, by=0.01)
# #x_points=x_points[-c(1,length(x_points))]

par1=c(1.5, 0.5, 0.5)
par2=c(1.5,1, 0.5)
par3=c(1.5, 1.5, 0.5)
par4=c(1.5, 2, 0.5)
par5=c(1.5, 2.5, 0.5)


(g1=deb_fast(par=par1, x=x_points))
(g2=deb_fast(par=par2, x=x_points))
(g3=deb_fast(par=par3, x=x_points))
(g4=deb_fast(par=par4, x=x_points))
(g5=deb_fast(par=par5, x=x_points))


plot(x_points, g1,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, rho=0.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
) 

lines(x_points, g2, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$beta$"), 
       c(0.5,1,1.5,2, 2.5), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)

#################################
### Vatiando rho
# x11()
# par(mfrow=c(1,2)) 

# x_points = seq(0,1, by=0.01)
# #x_points=x_points[-c(1,length(x_points))]

par1=c(1.5, 2.5, -0.95)
par2=c(1.5,2.5, -0.5)
par3=c(1.5, 2.5,0)
par4=c(1.5, 2.5, 0.5)
par5=c(1.5, 2.5, 0.95)


(g1=deb_fast(par=par1, x=x_points))
(g2=deb_fast(par=par2, x=x_points))
(g3=deb_fast(par=par3, x=x_points))
(g4=deb_fast(par=par4, x=x_points))
(g5=deb_fast(par=par5, x=x_points))


plot(x_points, g1,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, beta=2.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#,aumentar o tamanho do titulo
     #cex.lab = 1.5, cex.axis = 1.5
) 

lines(x_points, g2, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.95,-0.5,0,0.5, 0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)

#############################################
## investigando simetria e bimodalidade

x11()
par(mfrow=c(2,4)) 

x_points = seq(0,1, by=0.01)
#x_points=x_points[-c(1,length(x_points))]

par1=c(1.5, 0.5,-0.5)
par2=c(1.5, 0.5,-0.65)
par3=c(1.5, 0.5,-0.75)
par4=c(1.5, 0.5,-0.85)
par5=c(1.5, 0.5,-0.95)


(g11=deb_fast(par=par1, x=x_points))
(g21=deb_fast(par=par2, x=x_points))
(g31=deb_fast(par=par3, x=x_points))
(g41=deb_fast(par=par4, x=x_points))
(g51=deb_fast(par=par5, x=x_points))


plot(x_points, g11,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, beta=0.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g21, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g31, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g41, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g51, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)


#################################
### Variando 
par1=c(1.5, 1.5, -0.5)
par2=c(1.5,1.5, -0.65)
par3=c(1.5, 1.5, -0.75)
par4=c(1.5, 1.5, -0.85)
par5=c(1.5, 1.5, -0.95)


(g12=deb_fast(par=par1, x_points))
(g22=deb_fast(par=par2, x_points))
(g32=deb_fast(par=par3, x_points))
(g42=deb_fast(par=par4, x_points))
(g52=deb_fast(par=par5, x_points))


plot(x_points, g12,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, beta=1.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo)
)

lines(x_points, g22, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g32, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g42, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g52, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)


####################################


par1=c(1.5, 3.5, -0.5)
par2=c(1.5,3.5, -0.65)
par3=c(1.5, 3.5,-.75)
par4=c(1.5, 3.5, -0.85)
par5=c(1.5, 3.5, -0.95)


(g13=deb_fast(par=par1, x_points))
(g23=deb_fast(par=par2, x_points))
(g33=deb_fast(par=par3, x_points))
(g43=deb_fast(par=par4, x_points))
(g53=deb_fast(par=par5, x_points))


plot(x_points, g13,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, beta=2.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#,aumentar o tamanho do titulo)
)

lines(x_points, g23, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g33, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g43, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g53, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)





#################################
### Variando 

par1=c(1.5, 3.5, -0.5)
par2=c(1.5,3.5, -0.65)
par3=c(1.5, 3.5,-0.75)
par4=c(1.5, 3.5, -0.85)
par5=c(1.5, 3.5, -0.95)


(g14=deb_fast(par=par1, x_points))
(g24=deb_fast(par=par2, x_points))
(g34=deb_fast(par=par3, x_points))
(g44=deb_fast(par=par4, x_points))
(g54=deb_fast(par=par5, x_points))


plot(x_points, g14,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, beta=3.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#,aumentar o tamanho do titulo)
)

lines(x_points, g24, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g34, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g44, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g54, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)


#############################################
## investigando simetria e bimodalidade

# x11()
# par(mfrow=c(1,4)) 

#x_points = seq(0,1, by=0.01)
#x_points=x_points[-c(1,length(x_points))]

par1=c(0.5, 2.5,-0.5)
par2=c(0.5, 2.5,-0.65)
par3=c(0.5, 2.5,-0.75)
par4=c(0.5, 2.5,-0.85)
par5=c(0.5, 2.5,-0.95)


(g15=deb_fast(par=par1, x_points))
(g25=deb_fast(par=par2, x_points))
(g35=deb_fast(par=par3, x_points))
(g45=deb_fast(par=par4, x_points))
(g55=deb_fast(par=par5, x_points))


plot(x_points, g15,
     ylim = c(0,8),
     main=TeX("$alpha=0.5, beta=2.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g25, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g35, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g45, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g55, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)


#################################
### Variando 
par1=c(1.0, 2.5, -0.5)
par2=c(1.0,2.5, -0.65)
par3=c(1.0, 2.5, -0.75)
par4=c(1.0, 2.5, -0.85)
par5=c(1.0, 2.5, -0.95)


(g16=deb_fast(par=par1, x_points))
(g26=deb_fast(par=par2, x_points))
(g36=deb_fast(par=par3, x_points))
(g46=deb_fast(par=par4, x_points))
(g56=deb_fast(par=par5, x_points))


plot(x_points, g16,
     ylim = c(0,8),
     main=TeX("$alpha=1.0, beta=2.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo)
)

lines(x_points, g26, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g36, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g46, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g56, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)


####################################


par1=c(1.5, 2.5, -0.5)
par2=c(1.5,2.5, -0.65)
par3=c(1.5, 2.5,-0.75)
par4=c(1.5, 2.5, -0.85)
par5=c(1.5, 2.5, -0.95)


(g17=deb_fast(par=par1, x_points))
(g27=deb_fast(par=par2, x_points))
(g37=deb_fast(par=par3, x_points))
(g47=deb_fast(par=par4, x_points))
(g57=deb_fast(par=par5, x_points))


plot(x_points, g17,
     ylim = c(0,8),
     main=TeX("$alpha=1.5, beta=2.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#,aumentar o tamanho do titulo)
)

lines(x_points, g27, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g37, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g47, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g57, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95),
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)





#################################
### Variando 

par1=c(1.75, 2.5, -0.5)
par2=c(1.75,2.5, -0.65)
par3=c(1.75, 2.5,-0.75)
par4=c(1.75, 2.5, -0.85)
par5=c(1.75, 2.5, -0.95)


(g18=deb_fast(par=par1, x_points))
(g28=deb_fast(par=par2, x_points))
(g38=deb_fast(par=par3, x_points))
(g48=deb_fast(par=par4, x_points))
(g58=deb_fast(par=par5, x_points))


plot(x_points, g18,
     ylim = c(0,8),
     main=TeX("$alpha=1.75, beta=2.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#,aumentar o tamanho do titulo)
)

lines(x_points, g28, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g38, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g48, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g58, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$rho$"), 
       c(-0.5,-0.65,-0.75,-0.85, -0.95), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)

#############################################
## investigando simetria 

### plots com alpha=beta

x11()
par(mfrow=c(2,2)) 

x_points = seq(0,1, by=0.01)
#x_points=x_points[-c(1,length(x_points))]

par1=c(0.5, 0.5,0.5)
par2=c(0.75, .75,0.5)
par3=c(1.0, 1.0,0.5)
par4=c(1.5, 1.5,0.5)
par5=c(1.75, 1.75,0.5)


(g1a=deb_fast(par=par1, x_points))
(g2a=deb_fast(par=par2, x_points))
(g3a=deb_fast(par=par3, x_points))
(g4a=deb_fast(par=par4, x_points))
(g5a=deb_fast(par=par5, x_points))


plot(x_points, g1a,
     ylim = c(0,8),
     main=TeX("$rho=0.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g2a, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3a, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4a, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5a, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$alpha=beta$"), 
       c(0.5,0.75,1.0,1.5, 1.75), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)



par1=c(0.5, 0.5,0.75)
par2=c(0.75, .75,0.75)
par3=c(1.0, 1.0,0.75)
par4=c(1.5, 1.5,0.75)
par5=c(1.75, 1.75,0.75)


(g1b=deb_fast(par=par1, x_points))
(g2b=deb_fast(par=par2, x_points))
(g3b=deb_fast(par=par3, x_points))
(g4b=deb_fast(par=par4, x_points))
(g5b=deb_fast(par=par5, x_points))


plot(x_points, g1b,
     ylim = c(0,8),
     main=TeX("$rho=0.75$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g2b, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3b, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4b, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5b, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$alpha=beta$"), 
       c(0.5,0.75,1.0,1.5, 1.75), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.5#controla o tamanho da legenda) 
)

par1=c(0.5, 0.5,-0.5)
par2=c(0.75, .75,-0.5)
par3=c(1.0, 1.0,-0.5)
par4=c(1.5, 1.5,-0.5)
par5=c(1.75, 1.75,-0.5)


(g1c=deb_fast(par=par1, x_points))
(g2c=deb_fast(par=par2, x_points))
(g3c=deb_fast(par=par3, x_points))
(g4c=deb_fast(par=par4, x_points))
(g5c=deb_fast(par=par5, x_points))


plot(x_points, g1c,
     ylim = c(0,8),
     main=TeX("$rho=-0.5$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g2c, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3c, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4c, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5c, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$alpha=beta$"), 
       c(0.5,0.75,1.0,1.5, 1.75), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.2#controla o tamanho da legenda) 
)



par1=c(0.5, 0.5,-0.75)
par2=c(0.75, .75,-0.75)
par3=c(1.0, 1.0,-0.75)
par4=c(1.5, 1.5,-0.75)
par5=c(1.75, 1.75,-0.75)


(g1d=deb_fast(par=par1, x_points))
(g2d=deb_fast(par=par2, x_points))
(g3d=deb_fast(par=par3, x_points))
(g4d=deb_fast(par=par4, x_points))
(g5d=deb_fast(par=par5, x_points))


plot(x_points, g1d,
     ylim = c(0,8),
     main=TeX("$rho=-0.75$"),
     ylab = TeX("$f_Z(z)$"),xlab ="z" ,type="l",
     pch=15, lty=1, lwd=2.0,
     cex.main = 1.5#aumentar o tamanho do titulo
)

lines(x_points, g2d, type="l", 
      pch=16, 
      lty=2, 
      lwd=2.0,
      col="red") 
lines(x_points, g3d, type="l", 
      pch=17, 
      lty=3, 
      lwd=2.0,
      col="blue") 
lines(x_points, g4d, type="l", 
      pch=18, 
      lty=4, 
      lwd=2.0,
      col="gray") 
lines(x_points, g5d, type="l", 
      pch=19, 
      lty=5, 
      lwd=2.0,
      col="green") 


legend("topright", #inset=.05,
       y=0.5,
       title=TeX("$alpha=beta$"), 
       c(0.5,0.75,1.0,1.5, 1.75), 
       lty=1:5, col=c("black",  "red","blue", "gray", "green"),#, "magenta","purple","orange","yellow" ),
       cex=1.2#controla o tamanho da legenda) 
)







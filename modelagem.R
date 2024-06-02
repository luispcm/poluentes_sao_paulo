library(tidyverse)
library(MASS)
library(lubridate)
library(invgamma)
library(mvtnorm)
library(mice)
library(coda)
library(tsibble)
library(fpp3)
library(forecast)

set.seed(9)

base_sp = base %>% 
  mutate(ano = year(data),
         mes = month(data),
         sem = week(data)) %>% 
  group_by(ano,sem, região, lat, long) %>% 
  summarise(temp_media = mean(temp_media, na.rm = T),
            co_medio = mean(co_medio,na.rm = T),
            mp2.5_medio = mean(mp2.5_medio,na.rm = T)) %>% 
  mutate(across(everything(), ~replace(., is.nan(.), NA))) %>% 
  filter(ano == 2022,
         região == "Cid.Universitária-USP-Ipen"|
         região == "Congonhas"|
         região == "Ibirapuera"|
         região == "Itaim Paulista"|
         região == "Marg.Tietê-Pte Remédios"|
         região == "Mooca"|
         região == "Parque D.Pedro II"|
         região == "Perus"|
         região == "Pinheiros"|
         região == "Santana")



#definindo Y
y = c()

for (i in 1:length(unique(base_sp$sem))){
  
  y = rbind(y,base_sp$mp2.5_medio[base_sp$sem == i])
  
}

imputed_model = mice(y, seed = 9)
new.y = complete(imputed_model)
Y = new.y
Y = as.matrix(Y)
######### Sazonalidade ########
#estimando sazonalidade
J = 53
tempo <- 1:J

SIN <- COS <- matrix(nrow=J,ncol=6)
for (i in 1:6) {
  COS[,i] <- cos(2*pi*i*tempo/52)
  SIN[,i] <- sin(2*pi*i*tempo/52)
}

fit.mp <- lm(Y ~ COS+SIN)
summary(fit.mp)
par(mfrow=c(1,2),mar=c(4,4,2,0.5),cex=1.2)
plot(fit.mp$res,ylab="residuals",main="Residuals vs. time")
qqnorm(fit.mp$res)
qqline(fit.mp$res)


mediafit.mp <- cbind(rep(1,J), COS,SIN)%*%fit.mp$coef
plot(Y[,10],type="l",xlab="weeks",ylab="mp",xaxt="n")
lines(as.vector(mediafit.mp[,10]), col = 'red')

plot(Y[,4],type="l",xlab="weeks",ylab="mp",xaxt="n")
lines(as.vector(mediafit.mp[,4]), col = 'red')

#retirando sazonalidade
Y.final = Y - matrix(rep(mediafit.mp,n),J,n)
plot(Y.final[,8],type="l",xlab="weeks",ylab="mp",xaxt="n")
abline(h = 0, col = 'red')

#autocorrelação
acf(Y.final[8,],main="",ylab="ACF - mp")

Y = t(Y.final) #tem como desfazer, basta rodas a linha 48 e 49

####### gerando matrix de variáveis explicativas: ########
n = length(unique(base_sp$região))
lat = base_sp$lat
long = base_sp$long

X.mat = as.matrix(cbind(rep(1,n),
                  lat,
                  long))
X.mat = X.mat[1:10,]

#construindo matriz de distâncias:
matriz_longlat = cbind(X.mat[,2],X.mat[,3])

H = dist(matriz_longlat)
H <- as.matrix(H)


#####################################################################
#####################################################################
SOMA <- function(Rep,MAT,Y,mu){
  TOTAL <- 0
  for(k in 1:Rep){
    soma <- dmvnorm(Y[,k],mu,MAT,log=TRUE)
    TOTAL <- TOTAL + soma}
  return(TOTAL)
}
SOMA.2 = function(t,Y,R.inv,mu.c){
  total = 0
  for (i in 1:t){
    dif = Y[,i] - mu.c
    soma = (t(dif) %*% (R.inv) %*% (dif))/2
    total = total + soma
  }
  return(total)
}
### MCMC

#priori pra beta (intercepto, lat, long)

b = c(1,1,1)
B <- diag(100,3)

#priori pra sigma2

a0 <- 1
a1 <- 0.01
curve(dinvgamma(x,a0,a1),from = 0, to = 100)

#priori pra phi

d1 = 1
d0 = d1*median(H)
curve(dgamma(x,d0,d1),from = 0, to = 10)

# Tamanho da cadeia
M <- 20000

# Vetores das cadeias

beta.cadeia = matrix(nrow = M,
                     ncol = 3)

sigma2.cadeia = rep(NA,M)

phi.cadeia = rep(NA,M)

# Valores iniciais

beta.cadeia[1,] = c(0,0,-0.2)
sigma2.cadeia[1] = 34
phi.cadeia[1] = 0.175

## Contas previas
B.inv = solve(B)

contador = 0
R.c = exp(-H/phi.cadeia[1])
R.inv = solve(R.c)
time = length(unique(base_sp$sem))
Y.soma = rowSums(Y)

## Geração da amostra
for (i in 2:M){
  
  # Bloco BETA (via Gibbs)
  
  cov.inv <- 1/sigma2.cadeia[i-1]*R.inv
  
  produto <- t(X.mat)%*%cov.inv
  
  cov.beta = solve(time*produto%*%X.mat+B.inv)
  media.beta = cov.beta%*%(produto%*%Y.soma+B.inv%*%b)
  
  beta.cadeia[i,] = mvrnorm(1, media.beta, cov.beta)

  mu.c <- X.mat%*%beta.cadeia[i,]
  
  # Bloco SIGMA2 (via Gibbs)
  
  a0.star = a0 + n*time/2
  a1.star = SOMA.2(time,Y,R.inv,mu.c)  + a1
  
  sigma2.cadeia[i] = rinvgamma(1, a0.star, a1.star)
  
  vero.c = SOMA(time,sigma2.cadeia[i]*R.c,Y,mu.c)
  
  # Bloco PHI (via Metropoles Hastings)]
  # p significa "proposto" e c significa "cadeia"
  
  phi.p = exp(rnorm(1,log(phi.cadeia[i-1]),0.4))
  R.p = exp(-H/phi.p)
  vero.p = SOMA(time,sigma2.cadeia[i]*R.p,Y,mu.c)
  
  num = vero.p + dgamma(phi.p, d0,d1, log = T) + log(phi.p)
  den = vero.c + dgamma(phi.cadeia[i-1], d0,d1, log = T) + log(phi.cadeia[i-1])
  
  
  if (log(runif(1)) < (num - den)){
    phi.cadeia[i] <- phi.p
    contador <- contador + 1
    R.c <- R.p
    R.inv = solve(R.p)
  } else {
    phi.cadeia[i] <- phi.cadeia[i-1]
  }
}


t.a = contador*100/M
t.a


ini = 9

par(mfrow=c(3,2),mar=c(4,4,0.5,0.5))
plot(beta.cadeia[ini:M,1], type='l',ylab=expression(beta[0]))
plot(beta.cadeia[ini:M,2], type='l',ylab=expression(beta[lat]))
plot(beta.cadeia[ini:M,3], type='l',ylab=expression(beta[long]))
plot(sigma2.cadeia[ini:M], type='l',ylab=expression(sigma^2))
plot(phi.cadeia[ini:M], type='l',ylab=expression(phi))


######## Raftery Lewis Test #######

dado_exp = data.frame(sigma2.cadeia,beta.cadeia,phi.cadeia)
raftery.diag(dado_exp)

##### Amostra Efetiva ######
thinned_chain_exp = dado_exp[-c(seq(1,1000,1)),]
thin_factor <- 10
thinned_chain_exp <- thinned_chain_exp[c(seq(1, nrow(thinned_chain_exp), by = thin_factor)),]

par(mfrow=c(3,2),mar=c(4,4,0.5,0.5))
plot(thinned_chain_exp$X1, type='l',ylab=expression(beta[0]))
plot(thinned_chain_exp$X2, type='l',ylab=expression(beta[lat]))
plot(thinned_chain_exp$X3, type='l',ylab=expression(beta[long]))
plot(thinned_chain_exp$sigma2.cadeia, type='l',ylab=expression(sigma^2))
plot(thinned_chain_exp$phi.cadeia, type='l',ylab=expression(phi))

thinned_chain_exp = as.mcmc(thinned_chain_exp)
ic_exp = HPDinterval(thinned_chain_exp, prob = 0.95)
ic_exp = data.frame(ic_exp)
thinned_chain_exp = data.frame(thinned_chain_exp)

###### MODELO CAUCHY ############

# Vetores das cadeias

beta.cadeia = matrix(nrow = M,
                     ncol = 3)

sigma2.cadeia = rep(NA,M)

phi.cadeia = rep(NA,M)

# Valores iniciais

beta.cadeia[1,] = c(0,0,-0.2)
sigma2.cadeia[1] = 34
phi.cadeia[1] = 0.175

## Contas previas
B.inv = solve(B)

contador = 0
R.c = (1 + H/phi.cadeia[1])^(-1)
R.inv = solve(R.c)
time = length(unique(base_sp$sem))
Y.soma = rowSums(Y)

## Geração da amostra
for (i in 2:M){
  
  # Bloco BETA (via Gibbs)
  
  cov.inv <- 1/sigma2.cadeia[i-1]*R.inv
  
  produto <- t(X.mat)%*%cov.inv
  
  cov.beta = solve(time*produto%*%X.mat+B.inv)
  media.beta = cov.beta%*%(produto%*%Y.soma+B.inv%*%b)
  
  beta.cadeia[i,] = mvrnorm(1, media.beta, cov.beta)
  
  mu.c <- X.mat%*%beta.cadeia[i,]
  
  # Bloco SIGMA2 (via Gibbs)
  
  a0.star = a0 + n*time/2
  a1.star = SOMA.2(time,Y,R.inv,mu.c)  + a1
  
  sigma2.cadeia[i] = rinvgamma(1, a0.star, a1.star)
  
  vero.c = SOMA(time,sigma2.cadeia[i]*R.c,Y,mu.c)
  
  # Bloco PHI (via Metropoles Hastings)]
  # p significa "proposto" e c significa "cadeia"
  
  phi.p = exp(rnorm(1,log(phi.cadeia[i-1]),0.5))
  R.p = (1 + H/phi.p)^(-1)
  vero.p = SOMA(time,sigma2.cadeia[i]*R.p,Y,mu.c)
  
  num = vero.p + dgamma(phi.p, d0,d1, log = T) + log(phi.p)
  den = vero.c + dgamma(phi.cadeia[i-1], d0,d1, log = T) + log(phi.cadeia[i-1])
  
  
  if (log(runif(1)) < (num - den)){
    phi.cadeia[i] <- phi.p
    contador <- contador + 1
    R.c <- R.p
    R.inv = solve(R.p)
  } else {
    phi.cadeia[i] <- phi.cadeia[i-1]
  }
}


t.a = contador*100/M
t.a


ini = 500

par(mfrow=c(3,2),mar=c(4,4,0.5,0.5))
plot(beta.cadeia[ini:M,1], type='l',ylab=expression(beta[0]))
plot(beta.cadeia[ini:M,2], type='l',ylab=expression(beta[lat]))
plot(beta.cadeia[ini:M,3], type='l',ylab=expression(beta[long]))
plot(sigma2.cadeia[ini:M], type='l',ylab=expression(sigma^2))
plot(phi.cadeia[ini:M], type='l',ylab=expression(phi))


######## Raftery Lewis Test #######

dado_cauchy = data.frame(sigma2.cadeia,beta.cadeia,phi.cadeia)
raftery.diag(dado_cauchy)

##### Amostra Efetiva ######
thinned_chain_cauchy = dado_cauchy[-c(seq(1,1000,1)),]
thin_factor <- 10
thinned_chain_cauchy <- thinned_chain_cauchy[c(seq(1, nrow(thinned_chain_cauchy), by = thin_factor)),]

par(mfrow=c(3,2),mar=c(4,4,0.5,0.5))
plot(thinned_chain_cauchy$X1,type='l',ylab=expression(beta[0]))
plot(thinned_chain_cauchy$X2, type='l',ylab=expression(beta[lat]))
plot(thinned_chain_cauchy$X3, type='l',ylab=expression(beta[long]))
plot(thinned_chain_cauchy$sigma2.cadeia, type='l',ylab=expression(sigma^2))
plot(thinned_chain_cauchy$phi.cadeia, type='l',ylab=expression(phi))

thinned_chain_cauchy = as.mcmc(thinned_chain_cauchy)
ic_cauchy = HPDinterval(thinned_chain_cauchy, prob = 0.95)
ic_cauchy = data.frame(ic_cauchy)
thinned_chain_cauchy = data.frame(thinned_chain_cauchy)

#### Mediana e IC's ####

medianas = data.frame(
parametro = c("median_sigma2",
              "median_beta0",
              "median_beta1",
              "median_beta2",
              "median_phi"),
medianas_exp = c(median(thinned_chain_exp$sigma2.cadeia),
             median(thinned_chain_exp$X1),
             median(thinned_chain_exp$X2),
             median(thinned_chain_exp$X3),
             median(thinned_chain_exp$phi.cadeia)),
medianas_cauchy = c(median(thinned_chain_cauchy$sigma2.cadeia),
                              median(thinned_chain_cauchy$X1),
                              median(thinned_chain_cauchy$X2),
                              median(thinned_chain_cauchy$X3),
                              median(thinned_chain_cauchy$phi.cadeia)
))

######## Barra de erro ##########

ggplot(medianas,
       aes(x = parametro)) +
  geom_point(aes(x = parametro,
                 y = medianas_exp,
                 colour = "blue"),
             position = position_nudge(x = -0.1),
             size = 12) + 
  geom_point(aes(x = parametro,
                 y = medianas_cauchy,
                 colour = "red"),
             position = position_nudge(x = 0.1),
             size = 12) + 
  geom_errorbar(aes(ymin = ic_exp$lower, ymax = ic_exp$upper),
                width = 0.3,
                position = position_nudge(x = -0.1),
                size = 3) +
  geom_errorbar(aes(ymin = ic_cauchy$lower, ymax = ic_cauchy$upper),
                width = 0.3,
                position = position_nudge(x = 0.1),
                size = 3) +
  labs(x = "Parâmetros",
       y = "Mediana e I.C. de 95%",
       color = "Modelo") +
  theme_minimal() +
  scale_x_discrete(labels = c(expression(beta[0]),
                              expression(beta[1]),
                              expression(beta[2]),
                              expression(phi),
                              expression(sigma^2))) +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 26),
        legend.text=element_text(size=30),
        legend.title = element_text(size = 30),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c("Exponencial",
                                   "Cauchy")) +
  scale_y_continuous(
    limits = c(-21, 30),
    breaks = seq(-21, 30, by = 5)) +
  geom_hline(yintercept = 0,
             size = 1.5)
  
########## Gráfico dos Modelos #######
n1 = length(thinned_chain_exp$phi.cadeia)
rho = 0

H.simulado = seq(0, 0.5, length.out = 1000)

for (i in 1:n1){
  
  matriz = exp(-H.simulado/thinned_chain_exp$phi.cadeia[i])
  rho = rho + matriz
  
}

rho_medio_exp = (rho)/n1


rho = array(dim = c(1000,1,n1))
for (i in 1:n1){
      
      rho[,,i] = exp(-H.simulado/thinned_chain_exp$phi.cadeia[i])
}


rho_2.5_exp = apply(rho, MARGIN = c(1, 2), FUN = function(matri) quantile(matri, probs = 0.025))
rho_97.5_exp = apply(rho, MARGIN = c(1, 2), FUN = function(matri) quantile(matri, probs = 0.975))


################

n2 = length(thinned_chain_cauchy$phi.cadeia)
rho = 0

for (i in 1:n2){
  
  matriz = (1 + H.simulado/thinned_chain_cauchy$phi.cadeia[i])^(-1)
  rho = rho + matriz
  
}

rho_medio_cauchy = (rho)/n2


rho = array(dim = c(1000,1,n2))
for (i in 1:n2){
  
  rho[,,i] = (1 + H.simulado/thinned_chain_cauchy$phi.cadeia[i])^(-1)
}


rho_2.5_cauchy = apply(rho, MARGIN = c(1, 2), FUN = function(matri) quantile(matri, probs = 0.025))
rho_97.5_cauchy = apply(rho, MARGIN = c(1, 2), FUN = function(matri) quantile(matri, probs = 0.975))


######## Gráfico
ggplot(data = NULL,
       aes(x = H.simulado)) +
  geom_line(aes(y = rho_medio_exp,
            color = "rho_medio_exp"),
            size = 3) +
  geom_ribbon(aes(ymin = rho_2.5_exp,
                  ymax = rho_97.5_exp),
              fill = "blue",
              alpha = 0.5) +
  geom_line(aes(y = rho_medio_cauchy,
            color = "rho_medio_cauchy"),
            size = 3) +
  geom_ribbon(aes(ymin = rho_2.5_cauchy,
                  ymax = rho_97.5_cauchy),
              fill = "red",
              alpha = 0.5) +
  theme_minimal() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 45),
        legend.text=element_text(size= 30),
        legend.title = element_text(size = 30),
        legend.position = "bottom") +
  labs(x = "Distâncias",
       y = "Correlação Espacial a Posteriori",
       color = "Modelo") +
  scale_color_manual(values = c("rho_medio_exp" = "blue",
                                "rho_medio_cauchy" = "red"),
                     labels = c("Exponencial",
                                "Cauchy"))

##### AIC ######

AIC = c()
for (i in 1:n1){
  R = exp(-H/thinned_chain_exp$phi.cadeia[i])
  mu = X.mat%*%c(thinned_chain_exp$X1[i],
                 thinned_chain_exp$X2[i],
                 thinned_chain_exp$X3[i])
  vero[i] = SOMA(time,thinned_chain_exp$sigma2.cadeia[i]*R,Y,mu)
  AIC[i] = (-2) * vero[i] + 2 * 5
}

AIC_exp = mean(AIC)

AIC = c()
for (i in 1:n2){
  R = (1 + H/thinned_chain_cauchy$phi.cadeia[i])^(-1)
  mu = X.mat%*%c(thinned_chain_cauchy$X1[i],
                 thinned_chain_cauchy$X2[i],
                 thinned_chain_cauchy$X3[i])
  vero[i] = SOMA(time,thinned_chain_cauchy$sigma2.cadeia[i]*R,Y,mu)
  AIC[i] = (-2) * vero[i] + 2 * 5
}

AIC_cauchy = mean(AIC)

############## Retornanndo Y para a escala original #########

Y = new.y
Y = t(Y)


############ MAPE ##########
y_pred = array(NA, dim = c(10,53,n1))
mape = NULL
for (i in 1:n1){
  
  R = exp(-H/thinned_chain_exp$phi.cadeia[i])
  
  mu = X.mat%*%c(thinned_chain_exp$X1[i],
                 thinned_chain_exp$X2[i],
                 thinned_chain_exp$X3[i])
  
  y_pred[,,i] = t(mvrnorm(53, mu, thinned_chain_exp$sigma2.cadeia[i]*R))
  
  y_pred[,,i] = y_pred[,,i] + t(mediafit.mp)
    
  mape[i] = sum(abs((Y - y_pred[,,i])/Y))/530
}
mape_exp = mean(mape)*100

y_pred = array(NA, dim = c(10,53,n2))
mape = NULL
for (i in 1:n2){
  
  R = (1 + H/thinned_chain_cauchy$phi.cadeia[i])^(-1)
  
  mu = X.mat%*%c(thinned_chain_cauchy$X1[i],
                 thinned_chain_cauchy$X2[i],
                 thinned_chain_cauchy$X3[i])
  
  y_pred[,,i] = t(mvrnorm(53, mu, thinned_chain_exp$sigma2.cadeia[i]*R))
  
  y_pred[,,i] = y_pred[,,i] + t(mediafit.mp)
  
  mape[i] = sum(abs((Y - y_pred[,,i])/Y))/530
}

mape_cauchy = mean(mape)*100

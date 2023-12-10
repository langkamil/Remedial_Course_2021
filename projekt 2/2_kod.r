#install.packages('deSolve')

library(deSolve)

#LOGISTYCZNY
# N - wielkosc populacji w chwili t
# r - intensywnosc wzrostu populacji
# K - pojemnosc srodowiska


# Model deterministyczny  --------------------------------------------------


Logistyczny <-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # tempo zmiany
    dN <- r*N*(1-N/K)
    
    # wartosc tempa zmiany
    list(dN)
  }) 
}

#N=5 osobnikow, K=10 pojemnosc srodowiska, r - zmienia sie
times <- seq(0, 10, by = 0.01) #czas trwania
state <- c(N=5)

parameters_1 <- c(r=-1, K=10) 
out_1 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_1)
head(out_1)


parameters_2 <- c(r=-10, K=10) 
out_2 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_2)
head(out_2)


parameters_3 <- c(r=20, K=10) 
out_3 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_3)
head(out_3)


parameters_4 <- c(r=100, K=10) 
out_4 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_4)
head(out_4)

par(mfrow=c(2,2))
plot(out_1, ylab= "N", main="r = -1")
plot(out_2, ylab= "N", main="r = -10")
plot(out_3, ylab= "N", main="r = 20")
plot(out_4, ylab= "N", main="r = 100")

#N=5 osobnikow, r=10 wzrost, K - zmienia sie 

parameters_5 <- c(r=10, K=3) 
out_5 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_5)
head(out_5)


parameters_6 <- c(r=10, K=100) 
out_6 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_6)
head(out_6)


parameters_7 <- c(r=10, K=1000) 
out_7 <- ode(y = state, times = times, func = Logistyczny , parms = parameters_7)
head(out_7)

plot(out_5, ylab="N", main="K = 3")
plot(out_6, ylab="N", main="K = 100")
plot(out_7, ylab="N", main="K = 1000")

#r=10, K=100, N - zmienia sie 
parameters <- c(r=10, K=100) 
state <- c(N=5)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)
head(out)
plot(out, ylab="N", main="N = 5")

state <- c(N=50)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)
head(out)
plot(out, ylab="N", main="N = 50")

state <- c(N=100)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)
head(out)
plot(out, ylab="N", main="N = 100")


state <- c(N=200)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)
head(out)
plot(out, ylab="N", main="N = 200")


# _____________________________________________________________________
# Model stochastyczny -----------------------------------------------------

# przypadek (i)

intensywnosci1 <- function(i) {
  if (i <= 100  && i >= 0) return(c(i - i*i/100, i*i/100))
  else if(i > 100) return(c(0,  i*i/100))
}

# funckcja generujaca trajektorie dla przypadku (i)
traj1 <- function(n0, T) {
  n <- c(n0) # n zawiera liczbe osobnikow w danej chwili czasu t 
  t <- c(0) # czas 
  i <- 1
  while(t[i] < T && n[i] > 0) {
    
    # wyznaczamy intensywnosci
    intensywnosci <- intensywnosci1(n[i])
    
    # losujemy czas wykladniczy dla obu intensywnosci
    czas_lambda <- rexp(1, intensywnosci[1]) 
    czas_miu <- rexp(1, intensywnosci[2])
    
    # dokonujemy wyboru pomiedzy mniejsza z wartosci czasu 
    if(czas_lambda < czas_miu) {
      # jezeli czas dla intensywnosci lambda jest mniejszy to pojawia sie nowy osobnik
      t[i+1] <- t[i] + czas_lambda
      n[i+1] <- n[i] + 1
    } 
    else {
      # jezeli czas dla intensywnosci miu jest mniejszy to jeden osobnik z populacji umiera
      t[i+1] <- t[i] + czas_miu
      n[i+1] <- n[i] - 1
    }
    
    i <- i + 1
  }
  return(cbind(n,t))
}

# przypadek (ii)

intensywnosci2 <- function(i) {
  if (i >= 0) return(c(i, i*i/50))
}


# funckcja generujaca trajektorie dla przypadku (ii)
traj2 <- function(n0, T) {
  n <- c(n0) # n zawiera liczbe osobnikow w danej chwili czasu t 
  t <- c(0) # czas 
  i <- 1
  while(t[i] < T && n[i] > 0) {
    
    # wyznaczamy intensywnosci
    intensywnosci <- intensywnosci2(n[i])
    
    # losujemy czas wykladniczy dla obu intensywnosci
    czas_lambda <- rexp(1, intensywnosci[1]) 
    czas_miu <- rexp(1, intensywnosci[2])
    
    # dokonujemy wyboru pomiedzy mniejsza z wartosci czasu 
    if(czas_lambda < czas_miu) {
      # jezeli czas dla intensywnosci lambda jest mniejszy to pojawia sie nowy osobnik
      t[i+1] <- t[i] + czas_lambda
      n[i+1] <- n[i] + 1
    } 
    else {
      # jezeli czas dla intensywnosci miu jest mniejszy to jeden osobnik z populacji umiera
      t[i+1] <- t[i] + czas_miu
      n[i+1] <- n[i] - 1
    }
    
    i <- i + 1
  }
  return(cbind(n,t))
}


# porownanie wynikow deterministycznego i stochastycznego -------------

# (i)

x <- traj1(5,100)
times <- seq(0, 100, by = 0.01)
parameters <- c(r=1, K=50)
state <- c(N=5)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)

# wykres trajektorii procesu
plot(x[,2], x[,1], type='s', ylim=c(0,70), xlim=c(0,40), ylab="N(t)", xlab = "time")
lines(out[,1],out[,2],col=2)
legend(20, 20, legend=c("deterministyczny", "stochastyczny (i)"), col = c("red", "black"), lty=1:1, cex=0.8)

# porownanie wynikow deterministycznego i stochastycznego (ii)

z <- traj2(5,100)
times <- seq(0, 100, by = 0.01)
parameters <- c(r=1, K=50)
state <- c(N=5)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)

# wykres trajektorii procesu
plot(z[,2], z[,1], type='s', ylim=c(0,80), xlim=c(0,40),ylab="N(t)", xlab = "time")
lines(out[,1],out[,2],col=2)
legend(20, 20, legend=c("deterministyczny", "stochastyczny wersja (ii)"), col = c("red", "black"), lty=1:1, cex=0.8)



# Wyznaczenie wariancji oraz wartosci oczekiwanej w t = 10 dla (i) i (ii) ---------------------------


# (i)
wyniki_t10_i <- numeric(1000)
for (i in 1:1000) {
  x <- traj1(5,100)
  wyniki_t10_i[i] <-  x[floor(x[,2]) == 10, 1][1]
}
wartosc_oczekiwana_i <- mean(wyniki_t10_i)
wariancja_i <- var(wyniki_t10_i)


# (ii)
wyniki_t10_ii <- numeric(1000)
for (i in 1:1000) {
  z <- traj2(5,100)
  wyniki_t10_ii[i] <-  z[floor(z[,2]) == 10 ,1][1]
}
wartosc_oczekiwana_ii <- mean(wyniki_t10_ii)
wariancja_ii <- var(wyniki_t10_ii)


# funkcja wartosci oczekiwanej i wariancji --------------------------------

# (i)

sumy1 <- function(m, n0 , T, dt) {
  # podzial odcinka czasu (dyskretyzacja)
  dysk <- seq(0,T,by=dt)
  # suma trajektorii na podziale odcinka
  trajsuma <- numeric(length(dysk))
  # suma kwadratow trajektorii na podziale odcinka
  trajsuma2 <- numeric(length(dysk))
  
  # petla odpowiadajaca za trajektorie (bedzie ich m)
  for(i in 1:m){
    wynik <- traj1(n0, T)
    trajdysk <- sapply(dysk, function(x) mean(wynik[wynik[,2]<=x, 1]))
    trajsuma <- trajsuma + trajdysk
    trajsuma2 <- trajsuma2 + trajdysk^2
  }
  return(cbind(dysk,srednia=trajsuma/m,srednia2=trajsuma2/m))
}

# otrzymujemy dyskretny podzial czasu [0,T], srednia N(t) i srednia kwadratow (N^2(t))
T <- 100
dt <- 0.1
m <- 1000
n0 <- 5
wynik1 <- sumy1(m,n0,T,dt)

# wartosc oczekiwana
plot(wynik1[,1], wynik1[,2], type='l', xlab='t', ylab='E(N(T))')

# wariancja = (EN)^2 - (E(N^2))
plot(wynik1[,1], wynik1[,3]-wynik1[,2]^2,type='l',lwd=2,xlab='t',ylab='Var(N(T))')

# (ii)

sumy2 <- function(m, n0 , T, dt) {
  # podzial odcinka czasu (dyskretyzacja)
  dysk <- seq(0,T,by=dt)
  # suma trajektorii na podziale odcinka
  trajsuma <- numeric(length(dysk))
  # suma kwadratow trajektorii na podziale odcinka
  trajsuma2 <- numeric(length(dysk))
  
  # petla odpowiadajaca za trajektorie (bedzie ich m)
  for(i in 1:m){
    wynik <- traj2(n0, T)
    trajdysk <- sapply(dysk, function(x) mean(wynik[wynik[,2]<=x, 1]))
    trajsuma <- trajsuma + trajdysk
    trajsuma2 <- trajsuma2 + trajdysk^2
  }
  return(cbind(dysk,srednia=trajsuma/m,srednia2=trajsuma2/m))
}

# otrzymujemy dyskretny podzial czasu [0,T], srednia N(t) i srednia kwadratow (N^2(t))
T <- 100
dt <- 0.1
m <- 1000
n0 <- 5
wynik2 <- sumy2(m,n0,T,dt)

# wartosc oczekiwana
plot(wynik2[,1], wynik2[,2], type='l', xlab='t', ylab='E(N(T))')

# wariancja = (EN)^2 - (E(N^2))
plot(wynik2[,1], wynik2[,3]-wynik2[,2]^2,type='l',lwd=2,xlab='t',ylab='Var(N(T))')



# Porownanie E[N(t)] z rozwiazaniem deterministycznym ---------------------


# (i)

times <- seq(0, 100, by = 0.01)
parameters <- c(r=1, K=50)
state <- c(N=5)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)

plot(wynik1[,1], wynik1[,2], xlim=c(0,95), type='l', xlab='t', ylab='')
lines(out[,1],out[,2],col=2)
legend(60, 25, legend=c("deterministyczny: N(t)", "(i) stochastyczny: E[N(t)]"), col = c("red", "black"), lty=1:1, cex=0.8)

# (ii)

times <- seq(0, 100, by = 0.01)
parameters <- c(r=1, K=50)
state <- c(N=5)
out <- ode(y = state, times = times, func = Logistyczny , parms = parameters)

plot(wynik2[,1], wynik2[,2], xlim=c(0,40), type='l', xlab='t', ylab='')
lines(out[,1],out[,2],col=2)
legend(20, 25, legend=c("deterministyczny: N(t)", "(ii) stochastyczny: E[N(t)]"), col = c("red", "black"), lty=1:1, cex=0.8)

m <- 10000
pdf <- function(z1, z2) (1/(2*pi)*exp(-(z1^2 +z2^2)/2))
pdf_polar <- function(r) (1/(2*pi)*exp(-(r^2)/2)*r) # na podstawie rownania z rozwiazania analitycznego w zad. 2

# (a) P(0 < Z1 <= 1, 0 < Z2 <= 1) -------------------------------------------------

# MC ===================================

a1 <- 0
b1 <- 1 

a2 <- 0
b2 <- 1


MC_a <- sapply(1:100, function(x) sum(pdf(runif(m, a1, b1), runif(m, a2, b2))) * (b1 - a1)*(b2 - a2)/m)
MC_mean_a <- mean(MC_a)
MC_mean_a

# Riemann =================================== 

#podzial
x1 <- seq(a1, b1, by=1/m)
x2 <- seq(a2, b2, by=1/m)

R_a <-  sum(pdf(x1, x2))*(b1 - a1)^2/m
R_a

# (b) P(Z_1^2 + Z_2^2 < 1)  -----------------------------------------------

# MC ===================================

# biegunowe
r <- 1 * sqrt(runif(m))
phi <- runif(m) * 2 * pi

# kartezjanskie (kolo jednostkowe o srodku w (0,0))
s1 <- r * cos(phi)
s2 <- r * sin(phi)

MC_b <- sapply(1:100, function(i) sum(pdf(s1, s2)) * pi*1/m)
MC_mean_b <- mean(MC_b)
MC_mean_b


#Riemann ===================================

# na wspolrzednych biegunowych

# promien
r1 <- 0
r2 <- 1

# kat
phi1 <- 0
phi2 <- 2*pi

v1 <- seq(r1, r2, by=r2/m)

R_b <- sum(pdf_polar(v1))*(r2 - r1)*(phi2 - phi1)/m
R_b


# (c) P(Z1 > 0, Z2 > 0, Z1 + Z2  < 1) -----------------------------------------

# MC ===================================

a <- 0 
b <- 1
h=b-a

MC_c <- c()

for (j in 1:100) {
  
  MC_c_tmp <- c() 
  # generujemy liczby
  c1 <- a + (b - a) * runif(m) 
  c2 <- a + (b - a) * runif(m)
  
  for (i in 1:m) {
    
    if (c1[i] + c2[i] < 1) { # sprawdzamy czy sa w trojkacie
      
      MC_c_tmp = c(MC_c_tmp, pdf(c1[i], c2[i]))
      
    }
  }
  
  MC_c = c(MC_c, sum(MC_c_tmp)*h^2/m)
}

MC_mean_c <- mean(MC_c)
MC_mean_c



# Riemann ===================================


a <- 0
b <- 1 

k <- seq(a, b, 1/m)

R_c_tmp <- 0

for (j in k) {
  
  for(i in k) {
    
    if(j + i < 1) { # czy sa w trojkacie
      
      R_c_tmp <- R_c_tmp + pdf(j,i)
    }
  }
}

R_c <- R_c_tmp * (1/m)^2
R_c


# Eksperyment ------------------------------------------------------------

# Sprawdzmy dokladnosc wynikow dla wiekszego m

m_test <- 20000 

# (a)

# MC ===================================
MC_a_test <- sapply(1:100, function(x) sum(pdf(runif(m_test, a1, b1), runif(m_test, a2, b2))) * (b1 - a1)*(b2 - a2)/m_test)
MC_mean_a_test <- mean(MC_a_test)
MC_mean_a_test

# Riemann   ===================================
x1_test <- seq(a1, b1, by=1/m_test)
x2_test <- seq(a2, b2, by=1/m_test)

R_a_test <-  sum(pdf(x1_test, x2_test))*(b1 - a1)*(b2 - a2)/m_test
R_a_test


# (b)

# MC  ===================================

# biegunowe
r_test <- 1 * sqrt(runif(m_test))
phi_test <- runif(m_test) * 2 * pi

# kartezjanskie (kolo jednostkowe o srodku w (0,0))
s1_test <- r_test * cos(phi_test)
s2_test <- r_test * sin(phi_test)

MC_b_test <- sapply(1:100, function(i) sum(pdf(s1_test, s2_test)) * pi*1/m_test)
MC_mean_b_test <- mean(MC_b_test)
MC_mean_b_test


#Riemann  ===================================

# na wspolrzednych biegunowych

v1_test <- seq(r1, r2, by=r2/m_test)

R_b_test <- sum(pdf_polar(v1_test))*(r2 - r1)*(phi2 - phi1)/m_test
R_b_test


# (c)

# MC ===================================

MC_c_test <- c()

for (j in 1:100) {
  
  MC_c_tmp_test <- c()
  c1_test <- a + (b - a) * runif(m_test) 
  c2_test <- a + (b - a) * runif(m_test)
  
  for (i in 1:m_test) {
    
    if (c1_test[i] + c2_test[i] < 1) {
      
      MC_c_tmp_test = c(MC_c_tmp_test, pdf(c1_test[i], c2_test[i]))
      
    }
  }
  
  MC_c_test = c(MC_c_test, sum(MC_c_tmp_test)*h^2/m_test)
}

MC_mean_c_test <- mean(MC_c_test)
MC_mean_c_test



# Riemann ===================================


k_test <- seq(a, b, 1/m_test) # UWAGA: bardzo dlugo sie wykonuje

R_c_tmp_test <- 0

for (j in k_test) {
  
  for(i in k_test) {
    
    if(j + i < 1) {
      
      R_c_tmp_test <- R_c_tmp_test + pdf(j,i)
    }
  }
  print(j)
}

R_c_test <- R_c_tmp_test * (1/m_test)^2
R_c_test


summary <- t(data.frame( 'Metoda analityczna' = c(0.1165, 0.3935, NA,
                                                  0.1165, 0.3935, NA),
                         'Metoda Monte Carlo' = c(MC_mean_a, MC_mean_b, MC_mean_c,
                                                  MC_mean_a_test, MC_mean_b_test, MC_mean_c_test),
                         'Metoda Riemanna' = c(R_a, R_b, R_c,
                                               R_a_test, R_b_test, R_c_test),
                         row.names = c('(a) M = 10000', '(b) M = 10000', '(c) M = 10000',
                                       '(a) M = 20000', '(b) M = 20000', '(c) M = 20000')))
View(summary)

error <- data.frame('Blad wzgledny' = 
                      c(abs(0.1165 - MC_mean_a)*100/0.1165,
                        abs(0.3935 - MC_mean_b)*100/0.3935,
                        abs(NA - MC_mean_c)*100/NA,
                        abs(0.1165 - MC_mean_a_test)*100/0.1165,
                        abs(0.3935 - MC_mean_b_test)*100/0.3935,
                        abs(NA - MC_mean_c_test)*100/NA,
                        abs(0.1165 - R_a)*100/0.1165,
                        abs(0.3935 - R_b)*100/0.3935,
                        abs(NA - R_c)*100/NA,
                        abs(0.1165 - R_a_test)*100/0.1165,
                        abs(0.3935 - R_b_test)*100/0.3935,
                        abs(NA - R_c_test)*100/NA),
                    
                    row.names = c('(a) MC: M = 10000', '(b) MC: M = 10000', '(c) MC: M = 10000',
                                  '(a) MC: M = 20000', '(b) MC: M = 20000', '(c) MC: M = 20000',
                                  '(a) R: M = 10000', '(b) R: M = 10000', '(c) R: M = 10000',
                                  '(a) R: M = 20000', '(b) R: M = 20000', '(c) R: M = 20000'))
View(t(error))


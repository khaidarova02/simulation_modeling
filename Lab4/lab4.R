# Пункт 4
library(circular)
library(triangle)

# Загрузка данных
data <- read.delim("D:/ИмитМод/Lab4/PF_data.txt", sep=",")

# Инициализация параметров
theta_bar <- pi/3
L <- 2
s_bar <- 0.3
w_bar <- 0.1

N <- 1000
result <- c()
resultR <- c()

particles <- data.frame(theta = runif(N, min = 0, max = 2*pi),
                        R = runif(N, min = 0, max = 2*L))

# Функция для обновления весов частиц
update_weights <- function(z1, z2) {
  theta = particles$theta
  R = particles$R
  weights <- numeric(length(theta))

  for (i in 1:length(theta)) {
    z1_expected <- ifelse(theta[i] > theta_bar && theta[i] < 2*pi - theta_bar,
                          L - R[i]*cos(theta[i]),
                          L - R[i]*cos(theta_bar))
    z2_expected <- ifelse(theta[i] >= 0 && theta[i] <= pi, 1,-1)

    weight_z1 <- dtriangle(z1 - z1_expected, -w_bar, w_bar, 0)
    weight_z2 <- ifelse((is.na(z2) || z2 == z2_expected), 1, 0)

    weights[i] <- weight_z1 * weight_z2
  }
  weights <- weights / sum(weights)
  return(weights)
}

# Функция для ресэмплинга частиц
resample_particles <- function(particles, weights) {
  indices <- sample(1:nrow(particles), size = nrow(particles), replace = TRUE, prob = weights)
  particles <- particles[indices, ]
  return(particles)
}

# Функция для обновления theta
new_theta <- function(theta, R) {
  theta = particles$theta
  R = particles$R
  
  s <- runif(length(theta), -s_bar, s_bar)
  
  for (i in 1:length(theta)) {
    theta[i] <- ifelse(theta[i] > theta_bar && theta[i] < 2*pi - theta_bar,
                       (theta[i]+s[i]) %% (2*pi),
                       (theta[i]+ s[i]/(R[i]+s_bar)) %% (2*pi))
  }
  return(theta)
}

# Цикл обработки данных
for (k in 1:nrow(data)) {
  z1 <- data$distSensor[k]
  z2 <- data$halfPlaneSensor[k]
  
  # Обновление частиц
  particles$theta <- new_theta()
  
  # Обновление весов частиц
  weights <- update_weights(z1, z2)
  
  # Оценка текущего состояния
  result <- c(result, sum(particles$theta*weights))
  resultR <- c(resultR, sum(particles$R*weights))
  
  # Ресэмплинг частиц
  particles <- resample_particles(particles, weights)
}


# Визуализация результатов
plot((1:200), data$GroundTruth, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Theta")
lines((1:200), result, type = "l", col = "red", lwd = 2)
legend("topright", legend=c("Ground Truth", "Result Theta"), col=c("blue", "red"), lwd=2)
mean(resultR)

# Пункт 3.1

# Загрузка данных
data <- read.delim("D:/ИмитМод/Lab3/time_series.dat", sep=",")
data <- data$Y_n

# Инициализация параметров
theta <- 1
lambda <- 1
k <- 2
b1 <- 1
b2 <- 1

# Количество наблюдений
n <- length(data)

# Функция для генерации обратного гамма-распределения
rinvgamma <- function(alpha, beta) {
  return(1/rgamma(1, alpha, beta))
}

# Цикл Гиббса
n_iter <- 1000
samples <- matrix(NA, nrow = n_iter, ncol = 5)
for (iter in 1:n_iter) {
  
  # тут функции rgamma и rinvgamma были высчитаны - пример лекция стр. 136
  
  # Генерация theta из полного условного распределения
  theta <- rgamma(1, sum(data[1:k]), b1)
  
  # Генерация lambda из полного условного распределения
  lambda <- rgamma(1, sum(data[(k+1):n]), b2)
  
  # тут можно не брать min 
  
  # Генерация k из полного условного распределения с помощью алгоритма Метрополиса-Гастингса
  k_proposed <- sample(2:(n-1), 1)
  acceptance_prob <- min(1, sum(dpois(data[1:k_proposed], lambda, log = TRUE))
                         + sum(dpois(data[(k_proposed+1):n], theta, log = TRUE))
                         - sum(dpois(data[1:k], lambda, log = TRUE))
                         - sum(dpois(data[(k+1):n], theta, log = TRUE)))
  if (runif(1) < acceptance_prob) {
    k <- k_proposed
  }
  
  # Генерация b1 и b2 из обратного гамма-распределения
  b1 <- rinvgamma(0.5, theta)
  b2 <- rinvgamma(0.5, lambda)
  
  samples[iter, ] <- c(theta, lambda, k, b1, b2)
}

# Вывод результатов
colnames(samples) <- c("theta", "lambda", "k", "b1", "b2")
summary(samples)
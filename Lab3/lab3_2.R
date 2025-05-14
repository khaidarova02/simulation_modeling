# Пункт 3.2

# Загрузка данных
data1 <- read.delim("D:/ИмитМод/Lab3/test1.dat", sep=" ")
colnames(data1) <- c('doc','word','cnt_words')
data2 <- read.delim("D:/ИмитМод/Lab3/test2.dat", sep=" ")
colnames(data2) <- c('doc','word','cnt_words')

# Функция для collapsed Gibbs sampling
collapsed_gibbs_sampling <- function(data, T, alpha, beta, num_iter) {
  # Инициализация
  D <- max(data[,1])  # Количество документов
  V <- max(data[,2])  # Размер словаря
  N <- nrow(data)     # Общее количество слов
  
  # Инициализация счетчиков
  n_d_t <- matrix(0, nrow = D, ncol = T)  # Количество слов d, принадлежащих теме t
  n_t_w <- matrix(0, nrow = T, ncol = V)  # Количество слов из темы t, равных w
  n_d <- rep(0, D)  # Общее количество слов в документе d
  n_t <- rep(0, T)  # Общее количество слов в теме t
  
  z <- sample(1:T, N, replace = TRUE)  # Инициализация случайными темами
  
  for (i in 1:N) {
    d <- data[i, 1]
    w <- data[i, 2]
    t <- z[i]
    
    n_d_t[d, t] <- n_d_t[d, t] + data[i, 3]
    n_t_w[t, w] <- n_t_w[t, w] + data[i, 3]
    n_d[d] <- n_d[d] + data[i, 3]
    n_t[t] <- n_t[t] + data[i, 3]
  }
  
  # Гиббсовская выборка
  for (iter in 1:num_iter) {
    for (i in 1:N) {
      d <- data[i, 1]
      w <- data[i, 2]
      t <- z[i]
      
      # Убираем текущее слово из статистики
      n_d_t[d, t] <- n_d_t[d, t] - data[i, 3]
      n_t_w[t, w] <- n_t_w[t, w] - data[i, 3]
      n_d[d] <- n_d[d] - data[i, 3]
      n_t[t] <- n_t[t] - data[i, 3]
      
      # Вычисляем вероятности нового присвоения темы
      p_z <- rep(0, T)
      for (j in 1:T) {
        p_z[j] <- ((n_d_t[d, j] + alpha) / (n_d[d] + T * alpha)) * ((n_t_w[j, w] + beta) / (n_t[j] + V * beta))
      }
      p_z <- p_z / sum(p_z)
      
      # Выбираем новую тему
      new_t <- sample(1:T, 1, prob = p_z)
      
      # Обновляем статистику
      z[i] <- new_t
      n_d_t[d, new_t] <- n_d_t[d, new_t] + data[i, 3]
      n_t_w[new_t, w] <- n_t_w[new_t, w] + data[i, 3]
      n_d[d] <- n_d[d] + data[i, 3]
      n_t[new_t] <- n_t[new_t] + data[i, 3]
    }
  }
  
  return(z)
}

# Применение collapsed Gibbs sampling к данным test1.dat и test2.dat с заданными параметрами
z_test1 <- collapsed_gibbs_sampling(data1, T = 3, alpha = 1, beta = 1, num_iter = 100)
# второй код работает долго - Лукашенко его не попросит
# z_test2 <- collapsed_gibbs_sampling(data2, T = 20, alpha = 0.1, beta = 0.1, num_iter = 100)

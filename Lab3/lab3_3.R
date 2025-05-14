# Пункт 3.3

# здесь матрица итоговая как в лекциях нужна

# Инициализация параметров модели
beta <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3) # матрица коэффициентов взаимодействия
V <- 100                                                         # количество вершин графа
E <- matrix(c(1:V, 2:V, 1), ncol = 2)                            # определение ребер графа

# Функция для вычисления вероятности
compute_prob <- function(x, beta, E) {
  sum_exp <- 0
  for (edge in E) {
    i <- edge[1]
    j <- edge[2]
    sum_exp <- sum_exp + beta[x[i], x[j]]
  }
  return(exp(sum_exp))
}

# Инициализация случайного вектора X
X <- sample(1:3, V, replace = TRUE)

# Параметры алгоритма MCMC
n_iter <- 1000 # Количество итераций
accepted_samples <- matrix(NA, nrow = n_iter, ncol = V)

# Алгоритм Метрополиса-Гастингса
for (iter in 1:n_iter) {
  for (i in 1:V) {
    proposed_X <- X
    proposed_X[i] <- sample(setdiff(1:3, X[i]), 1)
    
    prob_ratio <- compute_prob(proposed_X, beta, E) / compute_prob(X, beta, E)
    acceptance_prob <- min(1, prob_ratio)
    
    if (!is.na(acceptance_prob) && runif(1) < acceptance_prob) {
      X <- proposed_X
    }
  }
  
  accepted_samples[iter, ] <- X
}

# Пример сгенерированных реализаций
print(accepted_samples[n_iter,])

# Сохранение результатов в файл
write.table(accepted_samples, file = "D:/ИмитМод/Lab3/accepted_samples.txt", sep = "\t", row.names = FALSE)

# Вывод сообщения об успешном сохранении
cat("Результаты сохранены в файл accepted_samples.txt\n")

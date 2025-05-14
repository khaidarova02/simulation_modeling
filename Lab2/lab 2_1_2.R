# Пункт 2.1.2
library(pracma)

# Инициализация переменных
center <- c(1,1,2)  # центр шара M
R <- 0.5              # радиус шара M

alpha <- 0.5          # некоторый положительный параметр для функции f(omega)

n <- 100              # количество экспериментов
count <- 1000         # количество векторов в эксперименте

# Функция микрофасетной BRDF
f <- function(theta) {
  return((alpha^2 * cos(theta)) / (pi * (cos(theta)^2 * (alpha^2 - 1) + 1)^2))
}

# Вспомогательная функция q1
q1 <- function(theta) {
  return(1/pi * cos(theta))
}

# Вспомогательная функция q2
q2 <- function(theta) {
  return(f(theta))
}

# Функция видимости V(omega)
V <- function(omega) {
  # Константы из уравнения,
  #   полученного в результате объединения уравнений прямой и шара
  a <- omega[1]^2 + omega[2]^2 + omega[3]^2
  b <- -2*(omega[1]*center[1] + omega[2]*center[2] + omega[3]*center[3])
  c <- center[1]^2 + center[2]^2 + center[3]^2 - R^2
  
  # Если дискриминант > 0, то пересекает шар,
  #   если = 0, то касается шара,
  #   иначе не пересекает и не касается шара
  discriminant <- b^2 - 4*a*c 
  if (discriminant >= 0) {
    return(0)
  }
  return(1)
}

# Функция для генерации случайных точек на полусфере
generate_points <- function(count_points, gen = "stand") {
  # Угол theta
  theta <- acos(1 - runif(count_points))            # равномерное распределение                
  if (gen == "q1")                              
  {
    theta <- acos(1-2*runif(count_points)) / 2      # распределение q1
  } else if (gen == "q2") {
    v <- runif(count_points)
    theta <- acos(sqrt( (1-v)/(v*(alpha^2-1)+1) ))  # распределение q2
  }
  
  # Угол phi
  phi <- runif(count_points, 0, 2 * pi)
  
  # Высчитывание новых координат
  x <- sin(theta) * cos(phi)
  y <- sin(theta) * sin(phi)
  z <- cos(theta)
  
  return(cbind(x, y, z))
}

# Метод Монте-Карло с использованием вспомогательного распределения q1
integral_q1 <- function() {
  points <- generate_points(count, "q1")
  # Считаем значения интеграла
  integrand <- apply(points, 1, function(omega) f(acos(omega[3])) * V(omega) * cos(acos(omega[3])) / q1(acos(omega[3])))
 
  return(mean(integrand))
}

# Метод Монте-Карло с использованием вспомогательного распределения q2
integral_q2 <- function() { 
  points <- generate_points(count, "q2")
  # Считаем значения интеграла
  integrand <- apply(points, 1, function(omega) f(acos(omega[3])) * V(omega) * cos(acos(omega[3])) / q2(acos(omega[3])))
  
  return(mean(integrand))
}

# Стандартный метод Монте-Карло
integral_standard <- function() {
  points <- generate_points(count)
  # Считаем значения интеграла
  integrand <- apply(points, 1, function(omega) f(acos(omega[3])) * V(omega) * cos(acos(omega[3])) * 2 * pi)
  
  return(mean(integrand))
}

# Пример оценки мат ожидания
integral_q1()
integral_q2()
integral_standard()

# Вычисляем дисперсии оценок
var_q1 <- var(replicate(n, integral_q1()))
var_q2 <- var(replicate(n, integral_q2()))
var_standard <- var(replicate(n, integral_standard()))

options(scipen = 1000)
cat("Дисперсия оценки интеграла с использованием вспомогательного распределения q1:", var_q1, "\n")
cat("Дисперсия оценки интеграла с использованием вспомогательного распределения q2:", var_q2, "\n")
cat("Дисперсия оценки интеграла с использованием стандартного метода Монте-Карло:", var_standard, "\n")

################################################################################

# Отрисовка распределений
Vdraw <- function(omega) {
  a <- omega[1]^2 + omega[2]^2 + omega[3]^2
  b <- -2*(omega[1]*center[1] + omega[2]*center[2] + omega[3]*center[3])
  c <- center[1]^2 + center[2]^2 + center[3]^2 - R^2
  
  discriminant <- b^2 - 4*a*c 
  if (discriminant > 0) {
    lines3d(x = c(0, omega[1]*((-b-sqrt(discriminant))/2*a)),
            y = c(0, omega[2]*((-b-sqrt(discriminant))/2*a)),
            z = c(0, omega[3]*((-b-sqrt(discriminant))/2*a)),
            col = "green", size = 3)
  } else if (discriminant == 0) {
    points3d(x = omega[1], y = omega[2], z = omega[3], col = "red", size = 3)
  } else {
    points3d(x = omega[1], y = omega[2], z = omega[3], col = "yellow", size = 3)
  }
}

draw <- function(points) {
  library(rgl)
  open3d()
  
  # Отрисовка координат
  axes3d(labels = c("X", "Y", "Z"))
  
  # Отображаем полусферу
  u <- seq(0, 2 * pi, length.out = 100)
  v <- seq(0, pi/2, length.out = 50)
  x <- 0 + 1 * outer(cos(u), sin(v))
  y <- 0 + 1 * outer(sin(u), sin(v))
  z <- 0 + 1 * outer(rep(1, length(u)), cos(v))
  surface3d(x, y, z, color = "blue", alpha = 0.3)
  
  # Отображаем шар
  u <- seq(0, 2 * pi, length.out = 100)
  v <- seq(0, pi, length.out = 50)
  x <- center[1] + R * outer(cos(u), sin(v))
  y <- center[2] + R * outer(sin(u), sin(v))
  z <- center[3] + R * outer(rep(1, length(u)), cos(v))
  surface3d(x, y, z, color = "red", alpha = 0.3)
  
  # Отображаем прямые
  w <- apply(points, 1, Vdraw)
}

draw(generate_points(n, "q1"))
draw(generate_points(n, "q2"))
draw(generate_points(n))

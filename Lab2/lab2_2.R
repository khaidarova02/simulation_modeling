# Пункт 2.2

# Инициализация переменных
N <- 10000        # количество фотонов
Np <- 0           # количество фотонов, прошедших через ткань
No <- 0           # количество отраженных фотонов

W_lim <- 10^(-4)  # пороговое значение веса 
m <- 10           # коэффициент нового веса
weight <- matrix(nrow=4, ncol=0)  # массив координат и поглощенного там веса

r <- 0.01         # радиус шара неоднородности
l <- 0.02         # глубина расположения шара неоднородности

paths <- list(list())  # список путей фотонов

data <- data.frame(
  name <- c("Эпидермис", "Папилл. дерма", "Поверхн. сосуд. сплетение",
            "Ретикул. дерма", "Глуб. сосуд. сплетение", "Неоднородность"),
  mu_a <- c(32, 23, 40, 23, 46, 51),
  mu_s <- c(165, 227, 246, 227, 253, 186),
  g <- c(0.72, 0.72, 0.72, 0.72, 0.72, 0.8),
  d <- c(0.01, 0.02, 0.02, 0.09, 0.06, 0)
)

# Вероятность поглащения
Pa <- function(layer) {
  return(data$mu_a[layer] / (data$mu_a[layer] + data$mu_s[layer]))
}

# Генерация угла theta при заданной плотности p(cos)
generate_theta <- function(layer, theta) {
  if (data$g[layer] != 0) {
    u <- (( 1 + data$g[layer]^2 - ((1 - data$g[layer]^2) / (1 - data$g[layer] + 2*data$g[layer]*theta))^2) / (2*data$g[layer]))
  } else {
    u <- 2 * theta - 1
  }
  return(acos(u))
}

# Генерация новых координат
new_coordinate <- function(layer, photon, V) {
  # Генерация пробега фотона
  dist <- rexp(1, data$mu_a[layer] + data$mu_s[layer])
  
  # Начальное погружение
  if (all(photon == c(0,0,0))) {
    return(list(c(0,0,dist), c(0,0,1)))
  }
  
  # Генерация углов отклонения траектории
  theta <- generate_theta(layer, runif(1))
  phi <- runif(1,0,2*pi)
  
  V <- V / sqrt(sum(V^2))
  if (abs(V[3]) != 1) {
    Vx = V[1]*cos(theta) + sin(theta) / sqrt(1-V[3]^2) * (V[1]*V[3]*cos(phi)-V[2]*sin(phi))
    Vy = V[2]*cos(theta) + sin(theta) / sqrt(1-V[3]^2) * (V[2]*V[3]*cos(phi)+V[1]*sin(phi))
    Vz = V[3]*cos(theta) - sin(theta)*cos(phi)*sqrt(1-V[3]^2)
  } else {
    Vx = sin(theta) * cos(phi)
    Vy = sin(theta) * sin(phi)
    Vz = V[3]*cos(theta)
  }
  
  # Высчитывание новых координат
  x = photon[1] + Vx * dist
  y = photon[2] + Vy * dist
  z = photon[3] + Vz * dist
  
  return(list(c(x,y,z), c(Vx,Vy,Vz)))
}

# Определение слоя
layer_definition <- function(dist){
  layer = 0
  if (dist < 0) {
    layer = -1
  } else if (dist < data$d[1]) {
    layer = 1
  } else if (dist < data$d[1] + data$d[2]) {
    layer = 2
  } else if (dist < data$d[1] + data$d[2] + data$d[1]) {
    layer = 3
  } else if (dist < data$d[1] + data$d[2] + data$d[3] + data$d[4]) {
    layer = 4
  } else if (dist < data$d[1] + data$d[2] + data$d[3] + data$d[4] + data$d[5]) {
    layer = 5
  } else {
    layer = 0
  }
  return(layer)
}

# Моделирование траекторий движения фотонов
for (i in 1:N) {
  layer <- 1               # начальный слой
  W <- 1                   # начальный статистический вес
  Wlost <- 0               # потерянный вес
  photon <- c(0,0,0)       # начальные координаты
  
  # Направляющий вектор
  vector <- c(0,0,1)
  
  # Запоминание координат для отображения пути фотона
  x <- c(0)
  y <- c(0)
  z <- c(0)
  while(TRUE) {
    # Перемещение фотона
    result <- new_coordinate(layer, photon, vector)
    photon <- result[[1]]
    vector <- result[[2]]
    
    # Определение текущего слоя
    layer <- layer_definition(photon[3])
    
    if (layer == -1 || layer == 0) {
      # Добавление координат и поглащенного веса
      weight <<- cbind(weight, c(photon[1], photon[2], photon[3], W))
      # Добавление пути текущего фотона
      paths[i] <- list(list(x,y,z))
      
      if (layer == -1) {
        # Подсчет отраженных фотонов
        No <- No + 1
      } else {
        # Подсчет фотонов, прошедших через ткань
        Np <- Np + 1
      }
      break;
    }
    
    # Проверка, что фотон находится в неоднородности
    dist <- sqrt(sum((photon - c(0,0,l))^2))
    if (dist <= r) {
      layer = 6
    }
    
    # Генерация дальнейшего действия над фотоном
    if (runif(1,0,1) < Pa(layer)) {
      # Добавление координат и поглащенного веса
      weight <<- cbind(weight, c(photon[1], photon[2], photon[3], W))
      # Добавление пути текущего фотона
      paths[i] <- list(list(x,y,z))
      break;
    }
    
    # Уменьшение статистического веса
    Wlost <- W*(data$mu_a[layer] / (data$mu_a[layer] + data$mu_s[layer]))
    W <- W - Wlost
    weight <<- cbind(weight, c(photon[1], photon[2], photon[3], Wlost))
    
    # Сравнение веса с пороговым
    if(W < W_lim) {
      if(runif(1,0,1) < 1/m) {
        # Добавление координат и поглащенного веса
        weight <<- cbind(weight, c(photon[1], photon[2], photon[3], W))
        # Добавление пути текущего фотона
        paths[i] <- list(list(x,y,z))
        break;
      }
      W = W * m
    }
    
    # Добавление координат в путь фотона
    x <- c(x, photon[1])
    y <- c(y, photon[2])
    z <- c(z, photon[3])
  }
}

cat("Коэффициент пропускания: ", Np/N, "\n")
cat("Коэффициент отражения: ", No/N, "\n")

################################################################################

# Поглащенная энергия в каждой слое биоткани
q <- N
Q <- matrix(0, nrow = 6, ncol = 2)
colnames(Q) <- c("V, см3", "Поглощенная энергия")

# Подсчет поглащенной энергии в слоях биоткани
z1 <- 0
for (i in 1:5) {
  V <- data$d[i]
  
  z2 <- z1 + data$d[i]
  X <- mean(weight[4, weight[3,] >= z1 & weight[3,] < z2
            & sqrt(sum((c(weight[1,],weight[2,],weight[3,]) - c(0,0,l))^2)) > r])
  if (is.nan(X)) {
    X <- 0
  }
  z1 <- z2
  
  Q[i,1] <- V
  Q[i,2] <- X*q / V
}

# Подсчет поглащенной энергии в неоднородности
V <- 4/3*pi*r^3
X <- mean(weight[4, sqrt(sum((c(weight[1,],weight[2,],weight[3,]) - c(0,0,l))^2)) <= r])
if (is.nan(X)) {
  X <- 0
}
Q[6,1] <- V
Q[6,2] <- X*q / V

# Вывод
Q

################################################################################

# Отрисовка распределения поглащенной энергии на оси xz
cnt_l <- 5
x1 <- -sum(data$d[1:cnt_l]) / 2
x2 <- sum(data$d[1:cnt_l]) / 2
z2 <- sum(data$d[1:cnt_l])
X <- weight[1, weight[1,] >= x1 & weight[1,] < x2 & weight[3,] >= 0 & weight[3,] < z2]
Z <- weight[3, weight[1,] >= x1 & weight[1,] < x2 & weight[3,] >= 0 & weight[3,] < z2]
W <- weight[4, weight[1,] >= x1 & weight[1,] < x2 & weight[3,] >= 0 & weight[3,] < z2] / max(weight[4, weight[1,] >= x1 & weight[1,] < x2 & weight[3,] >= 0 & weight[3,] < z2])
smoothScatter(X, Z, col = W, colramp=colorRampPalette(c("white", "red")),
              main = "Поглощенная энергия", xlab = "x, см", ylab = "z, см")

################################################################################

# Отрисовка прохождения фотонов через ткань
library(rgl)
open3d()

# Отрисовка координат
axes3d(labels = c("X", "Y", "Z"))

# Отрисовка слоев
x <- c(-0.1, 0.1, 0.1, -0.1)
y <- c(-0.1, -0.1, 0.1, 0.1)

z <- 0
quads3d(x=x, y=y, z=c(z,z,z,z), col="blue", alpha=0.5)
for (i in 1:(length(data$d)-1)) {
  z <- z + data$d[i]
  quads3d(x=x, y=y, z=c(z,z,z,z), col="blue", alpha=0.5)
}

# Отрисовка неоднородности
u <- seq(0, 2 * pi, length.out = 100)
v <- seq(0, pi, length.out = 50)
x <- 0 + r * outer(cos(u), sin(v))
y <- 0 + r * outer(sin(u), sin(v))
z <- l + r * outer(rep(1, length(u)), cos(v))
surface3d(x, y, z, color = "black", alpha = 0.3)

# Отрисовка фотонов
for (i in 1:N) {
  lines3d(x=paths[[i]][[1]], y=paths[[i]][[2]], z=paths[[i]][[3]], col="red")
}
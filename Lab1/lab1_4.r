# Пункт 1.4

# Инициализация переменных
R <- 15                 # Радиус блока в см
E <- 200                # Энергия выделяющаяся в одном акте деления в МэВ
p1 <- 0.7               # Вероятность захвата ядрами примесей
p2 <- 0.7               # Вероятность вылета 2-х нейтронов при делении ядер
mean_path <- 1          # Средний пробег нейтрона

n <- 100                # Количество итераций для метода Монте-Карло

# Вычиление расстояние нового нейтрона от центра блока
distance <- function(point) {
  return(sqrt(sum(point^2)))
}

# Генерация координат нового нейтрона
new_coordinate <- function(point) {
  # Генерация углов отклонения
  theta <- acos(1 - 2 * runif(1, 0, 1))
  phi <- runif(1, 0, 2 * pi)
  
  Vx = sin(theta) * cos(phi)
  Vy = sin(theta) * sin(phi)
  Vz = cos(theta)
  
  # Генерация пробега фотона
  s <- rexp(1, rate = 1/mean_path)
  
  # Высчитывание новых координат
  x = point[1] + Vx * s
  y = point[2] + Vy * s
  z = point[3] + Vz * s
  
  return(c(x,y,z))
}

# Рисование блока с нейтронами
draw_sphere <- function(neutrons) {
  library(rgl)
  open3d()
  
  # Рисуем шар
  u <- seq(0, 2 * pi, length.out = 100)
  v <- seq(0, pi, length.out = 50)
  x <- 0 + R * outer(cos(u), sin(v))
  y <- 0 + R * outer(sin(u), sin(v))
  z <- 0 + R * outer(rep(1, length(u)), cos(v))
  surface3d(x, y, z, color = "blue", alpha = 0.3)
  
  # Рисуем точки
  points3d(do.call(rbind, neutrons), col = "red", size = 5)
}

# Генерация появления нейтронов
neutron_cascade <- function(draw = FALSE)
{
  # Инициализация начальных параметров
  total_energy <- 0                   # Энергия, выделенная в блоке
  neutrons_list <- list(c(0, 0, 0))   # Список нейтронов
  
  # Набор нейтронов в блоке
  all_neutrons <- list(c(0, 0, 0))
  
  # Имитация развития нейтронного каскада
  while (length(neutrons_list) > 0) 
  {
    # Запоминаем координаты текущего нейтрона
    current_neutron <- neutrons_list[[1]]
    # Удаляем нейтрон из списка, тк он или поглатится, или поделится
    neutrons_list <- neutrons_list[-1]
    
    # Генерация события деления нейтрона
    if (runif(1) > p1) 
    {
      # Генерация количества появления вторичных нейтронов
      countn = 2 + rbinom(1, size = 1, prob = 1-p2)
      # Добавление выделевшейся энергии при акте деления
      total_energy <- total_energy + E
      
      # Генерация координатов для вторичных нейтронов
      for (i in 1:countn) {
        neutron_coordinate <- new_coordinate(current_neutron)
        if (distance(neutron_coordinate) <= R) {
          neutrons_list <- c(neutrons_list, list(neutron_coordinate))
          all_neutrons <- c(all_neutrons, list(neutron_coordinate))
          
          # Ограничение для случаев, когда множество нейтроном
          #   генерируются в блоке
          if (length(all_neutrons) == 100) {
            if (draw) {draw_sphere(all_neutrons)}
            return(total_energy)
          }
        }
      }
    }
  }
  if (draw) {draw_sphere(all_neutrons)}
  return(total_energy)
}

# Выполняем метод Монте-Карло для заданного количества итераций
energies <- replicate(n, neutron_cascade())

# Вычисляем среднее значение энергии
mean_energy <- mean(energies)
cat("Среднее количество энергии:", mean_energy, "МэВ\n")

# Пример сферы
neutron_cascade(TRUE)

Usage:

clang++ fluid.cpp -o fluid -std=c++20


## замена find на 4 ифа ( у нас всего 4 элемента)

С заменой 
Simulation completed in: 185.454 seconds

Без замены 
Simulation completed in: 678.596 seconds

# Перешел на g = 9.8, T = 1'000

## Добавление потоков для составления dirs
С заменой:
Simulation completed in: 7.97915 seconds

До этого:
Simulation completed in: 56.158 seconds

## Добавление потоков для velocity внутри цикла
С заменой:
Simulation completed in: 7.79796 seconds

Без замены: 
Simulation completed in: 7.97915 seconds

## Распаралеливание dfs

С заменой:
Simulation completed in: 2.22334 seconds

Без замены:
Simulation completed in: 7.79796 seconds

## Добавил тред пул 

Время особо не поменялось, но код стал чище :)
# Тестовое задание YADRO
специализация AME-L1 algorithms trainee 
<br/>
Выполнен второй вариант
Подробнее в папке [FFT](./FFT)
____

Была добавлена система сборки CMake 3.18 + Ninja  
И внесены изменения, чтобы собиралось с GCC     
В папке [scripts](./scripts) скрипты для простой сборки. В них уже указаны необходимые параметры, но можно при желании поменять
(Нужнен компилятор Clang или GCC, указаный в PATH)

``COMPILER_OPTION=clang|gcc``    
``*FLOAT_PREC=1|2``  
``PREALLOC_SIZE=(c++ uint32_t число) (по умолчанию 1<<15)``  
``*ARRAY_DFT_METHOD=ARRAY_DFT_ASYNC|ARRAY_DFT_THREADS``  
``NUM_DFT_THREADS=(c++ uint8_t число) (по умолчанию 4)``  

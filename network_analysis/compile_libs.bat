@gcc -o networks.o -c networks.c -L ./ -l_mt64
@ar rc lib_stat.a networks.o mt64.o
@gcc test.c -o test.exe -L ./ -l_stat

@echo Done



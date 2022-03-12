gcc -O3 -c rk.c ab.c DIYdefs.c
gcc -O3 DIYsimulate.c -o DIYsimulate.exe -l fftw3 rk.o ab.o DIYdefs.o
rm *.o

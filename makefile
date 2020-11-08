all: MD.exe
INCLUDE = ./

math.o: math.f90
	ifort -c math.f90 -lm -L/home/yongle/lib -lfftw3
main.o: main.f90
	ifort -c main.f90 system.f90 -lm 
MD.exe: main.o math.o system.o
	ifort -I/home/yongle/PIMD_1D/F90/800K_b -module . -o MD.exe main.o system.o math.o -lm -L/home/yongle/lib -lfftw3
system:
	ifort -c system.f90 -lm
clean:
	rm main.o math.o system.o MD.exe *.mod

file=ape_obs

all: build

build:
	gcc -c sizes.c -o sizes.o -Wall -Wextra
	ar rcs libsizes.a sizes.o
	gcc -c matops.c -o matops.o -Wall -Wextra -lm --openmp
	ar rcs libmatops.a matops.o
	gcc -c su3lib.c -o su3lib.o -Wall -Wextra -lm --openmp
	ar rcs libsu3lib.a su3lib.o
	gcc -c apelib.c -o apelib.o -Wall -Wextra -lm --openmp
	ar rcs libapelib.a apelib.o
	gcc $(file).c -I. -L. -lapelib -lsu3lib -lmatops -lsizes -llapacke -o $(file) -lm --openmp


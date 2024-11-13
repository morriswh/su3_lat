all: build

build:
	gcc -c sizes.c -o sizes.o -Wall -Wextra
	ar rcs libsizes.a sizes.o
	gcc -c su3lib.c -o su3lib.o -Wall -Wextra -lm --openmp
	ar rcs libsu3lib.a su3lib.o
	gcc main.c -I. -L. -lsu3lib -lsizes -o main -lm --openmp
	# gcc testing.c -I. -L. -lsu3lib -lsizes -o testing -lm --openmp


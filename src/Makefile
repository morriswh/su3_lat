all: build

build:
	gcc -c sizes.c -o sizes.o -Wall -Wextra
	ar rcs libsizes.a sizes.o
	gcc -c su3lib.c -o su3lib.o -Wall -Wextra -lm --openmp
	ar rcs libsu3lib.a su3lib.o
	gcc -c apelib.c -o apelib.o -Wall -Wextra -lm --openmp
	ar rcs libapelib.a apelib.o
	gcc main.c -I. -L. -lapelib -lsu3lib -lsizes -o main -lm --openmp
	# gcc testing.c -I. -L. -lsu3lib -lsizes -o testing -lm --openmp


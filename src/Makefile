# file=imPropNum
file=source_avg_num_less

all: build

build:
	gcc -c sizes.c -o sizes.o -Wall -Wextra
	ar rcs libsizes.a sizes.o
	gcc -c fileIO.c -o fileIO.o -Wall -Wextra
	ar rcs libfileIO.a fileIO.o
	gcc -c matops.c -o matops.o -Wall -Wextra -lm --openmp
	ar rcs libmatops.a matops.o
	gcc -c su3lib.c -o su3lib.o -Wall -Wextra -lm --openmp
	ar rcs libsu3lib.a su3lib.o
	gcc -c apelib.c -o apelib.o -Wall -Wextra -lm --openmp
	ar rcs libapelib.a apelib.o
	# gcc -c slat.c -I. -L. -lapelib -lsu3lib -lmatops -lsizes -lblas -llapack -o slat.o -Wall -Wextra
	# ar rcs libslat.a slat.o
	# gcc $(file).c -I. -L. -lapelib -lsu3lib -lmatops -lsizes -lblas -llapack -g -o $(file) -lm --openmp
	gcc $(file).c -I. -L. -lfileIO -lapelib -lsu3lib -lmatops -lsizes -llapack -g -o ../$(file) -lm --openmp
	# gcc $(file).c -I. -L. -lfileIO -lsu3lib -lsizes -g -o ../$(file) -lm --openmp

clean:
	rm *.o *.a


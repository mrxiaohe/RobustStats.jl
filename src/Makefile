all : akerd.so
 
akerd.so : akerd.h akerd.c
	gcc -O3 -fopenmp -fPIC -c akerd.c
	gcc -O3 -fopenmp -fPIC -shared -o libakerd.so akerd.o
	rm -f akerd.o
 
clean:
	rm -f akerd.so
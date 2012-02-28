
OBJFILES=util.o debug.o main.o drawfns.o triangulate.o circle_pack.o
DEPS =./include/* 

#COMPILE= g++ -Wall -I ./include -g -fno-inline -c $< -o $@ 
COMPILE= g++ -Wall -I ./include -g -c $< -o $@ 

all: $(OBJFILES)
	g++ -o draw $(OBJFILES) 

util.o: ./src/util.cpp ./include/util.h
	$(COMPILE)

debug.o: ./src/debug.cpp $(DEPS)
	$(COMPILE)


main.o: ./src/main.cpp $(DEPS)
	$(COMPILE)

drawfns.o: ./src/drawfns.cpp $(DEPS)
	$(COMPILE)

triangulate.o: ./src/triangulate.cpp $(DEPS)
	$(COMPILE)

circle_pack.o: ./src/circle_pack.c
	/usr/bin/gcc -g -O2 -s -fomit-frame-pointer -c $< -o $@

clean:
	rm -f $(OBJFILES)

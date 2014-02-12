CC = gcc
OPT = -W -Wall

OBJS = PMdef.o PMconjugation.o PMenlight.o PMextract.o PMdisplay.o \
       PMplanmap.o interface.o stats.o

all: $(OBJS) planarmap

PMdef.o: PMdef.c
	$(CC) $(OPT) -c PMdef.c -o PMdef.o

PMconjugation.o: PMconjugation.c
	$(CC) $(OPT) -c PMconjugation.c -o PMconjugation.o

PMenlight.o: PMenlight.c PMenlight.h
	$(CC) $(OPT) -c PMenlight.c -o PMenlight.o

PMextract.o: PMextract.c
	$(CC) $(OPT) -c PMextract.c -o PMextract.o

PMdisplay.o: PMdisplay.c PMdisplay.h
	$(CC) $(OPT) -c PMdisplay.c -o PMdisplay.o

PMplanmap.o: PMplanmap.c PMplanmap.h
	$(CC) $(OPT) -c PMplanmap.c -o PMplanmap.o

stats.o: stats.c stats.h
	$(CC) $(OPT) -c stats.c -o stats.o

planarmap: planarmap.c $(OBJS)
	$(CC) $(OPT) planarmap.c $(OBJS) -lm -o planarmap

interface.o: interface.c
	$(CC) $(OPT) -c interface.c -o interface.o

clean:
	/bin/rm -f $(OBJS)

VERSION=PlanarMap-v1-2b

distrib:
	tar zcvf $(VERSION).tgz *.c *.h Makefile README

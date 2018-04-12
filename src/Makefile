# Source, Executable, Includes, Library Defines
VERSION = 1.00
CC = gcc
INCL = kseq.h sabre.h
SRC = demulti_paired.c demulti_single.c sabre.c utils.c
OBJ = $(SRC:.c=.o)
DSRC=src

CFLAGS = -Wall -O2 -std=c99 -pedantic -DVERSION=$(VERSION)
LDFLAGS = -lz
GPROF = -pg
EXE = sabre

.PHONY: default

default: build
# a smarter way to have an if statement here instead of explicit grpof target
# have a look at gcc -M

%.o: %.c
	$(CC) -c $(CFLAGS) $(SRC)

demulti_single.o: kseq.h sabre.h
demulti_paired.o: kseq.h sabre.h
sabre.o: sabre.h

build: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS) 
	ln -sf $(DSRC)/$(EXE) ..

gprof: 
	$(CC) $(CFLAGS) $(GPROF) $(SRC) -o $(EXE).gprof $(LDFLAGS) 

clean:
	$(RM) $(OBJ) $(EXE) $(EXE).gprof core gmon.out

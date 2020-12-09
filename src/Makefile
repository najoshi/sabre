# Source, Executable, Includes, Library Defines
GIT_VERSION := "$(shell git describe --abbrev=7 --always --tags)"
CC = gcc
INCL = fastq.h sabre.h
#SRC = demulti_paired.c demulti_single.c sabre.c utils.c
SRC = sabre.c usage.c demultiplex.c utils.c fastq.c
OBJ = $(SRC:.c=.o)
DSRC=src

#CFLAGS = -Wall -O2 -std=c99 -pedantic -DVERSION=$(VERSION)
# need to quote GIT_VERSION so that the value gets passed as a string
CFLAGS = -Wall -O2 -std=gnu99 -pedantic -DVERSION=\"$(GIT_VERSION)\"
CFLAGSDEV = -Wall -O0 -g -std=gnu99 -DVERSION=\"$(GIT_VERSION)-dev\"

LDFLAGS = -lz -lpthread
GPROF = -pg
EXE = sabre

.PHONY: default

default: build
# a smarter way to have an if statement here instead of explicit grpof target
# have a look at gcc -M

%.o: %.c
	$(CC) -c $(CFLAGS) $(SRC)

usage.o: usage.h sabre.h fastq.h
utils.o: utils.h sabre.h fastq.h
demultiplex.o: demultiplex.h utils.h sabre.h fastq.h
sabre.o: sabre.h

build: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS)
	#ln -sf $(DSRC)/$(EXE) ../sabre

dev: $(OBJDEV)
	$(CC) $(CFLAGSDEV) $(SRC) -o $(EXE)-dev $(LDFLAGS)

metrics:
	$(CC) $(CFLAGSDEV) -o metrics metrics.c $(LDFLAGS)

gprof:
	$(CC) $(CFLAGS) $(GPROF) $(SRC) -o $(EXE).gprof $(LDFLAGS)

clean:
	$(RM) $(OBJ) $(EXE) core

clean-all:
	$(RM) $(OBJ) $(EXE) $(EXE)-dev $(EXE).gprof core gmon.out metrics

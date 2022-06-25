IDIR = include  -I/home/vorosadam/"Numerical Recipes in C"/nr_c304/code
ODIR = obj
BDIR = bin
SDIR = src
CC = g++
CFLAGS_R = -I$(IDIR) -lm
CFLAGS_C = -I$(IDIR)  

all:$(BDIR)/hatar

$(BDIR)/hatar: $(ODIR)/hatar.o
	$(CC) -o $@  $< $(CFLAGS_R)

$(ODIR)/%.o: $(SDIR)/%.C 
	$(CC) -o $@ -c $< $(CFLAGS_C)

.PHONY: clean

clean: 
	rm $(BDIR)/* $(ODIR)/*

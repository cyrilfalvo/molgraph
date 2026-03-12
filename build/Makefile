#
#  Makefile 
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
#


CC       = icpc
LD       = icpc

CFLAGS   = -Wall -O3
#CFLAGS   = -g
LIBS     = -lnauty

#
PROG     = ../molgraph.x
LISTEOBJ = \
        main.o

#
# compilation
#
%.o: ../src/%.cpp
	$(CC) -c $(CFLAGS) -o $@ $<

#
# dependances
#
all : $(PROG)
#
# edition de lien
#
$(PROG) : $(LISTEOBJ)
	$(LD) -o $(PROG) $(CFLAGS) $(LISTEOBJ) $(LIBS)
#
# effacement des fichiers intermediaires
#
clean :
	rm -f $(PROG) $(LISTEOBJ) core a.out



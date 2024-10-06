include ../make.inc
TSPACEHOME=http://aquarius.mp.es.osaka-u.ac.jp/~tspace/tspace_main/tsp07
WGET=wget
all: tsp98.f prmtsp.f tsp98.o
tsp98.f:
	$(WGET) $(TSPACEHOME)/$@
prmtsp.f:
	$(WGET) $(TSPACEHOME)/$@
clean:
	rm tsp98.o

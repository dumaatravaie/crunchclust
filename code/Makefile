CFLAGS = -O3
SRC=CrunchClust_Version43.cpp

all: clean crunchclust

crunchclust: $(SRC) 
	g++ $(CFLAGS) $^ -o $@

clean: 
	rm -f *.o crunchclust 


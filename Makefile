DEP = fm_random.o fuzzymeasuretools.o fmrandexample.cpp   

all: fmrandexample

fmrandexample:	$(DEP) 	
	g++ -O3 fmrandexample.cpp fuzzymeasuretools.o fm_random.o -o fmrandexample	

fm_random.o: 
	g++ -O3 -c fm_random.cpp -o fm_random.o

fuzzymeasuretools.o: 
	g++ -O3 -c fuzzymeasuretools.cpp -o fuzzymeasuretools.o

clean:
	rm -rf *.o
	rm fmrandexample



# Hvb makefile

hvb=hvb121212/
extra=lib/

FLAGS=-pedantic -Wall -g -I$(hvb) -I$(extra) -O3 
LIBS=-llapack -lblas -lm -L/usr/X11R6/lib -lX11 -lImlib2 

HVBS= $(hvb)matrix.o $(hvb)common.o $(hvb)easyx.o $(hvb)easyps.o $(hvb)graph.o $(hvb)parsing.o $(hvb)cmatrix.o $(hvb)easyim.o $(hvb)fitting.o $(hvb)optimize.o $(hvb)plotting.o $(hvb)manybody.o

EXTRA= $(extra)wf.o $(extra)complexplot.o

objects: $(HVBS) $(EXTRA) $(OBJECTS)

%.o: %.cc
	g++ -c $< $(FLAGS) -o $@

x%: x%.o objects
	g++ $(FLAGS) -o $@ $< $(HVBS) $(EXTRA) $(LIBS)

clean:
	ls x* | grep -v "\." | xargs rm -f
	rm *.o 

mrproper: 
	rm $(hvb)/*.o $(extra)/*.o


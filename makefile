F77 = mpif90 -i4 -real-size 32 -O4

FILES = maind.f gutsf.f gutsp.f misc.f boundary.f part_init.f initial.f
INCLUDE = incurv.h para.h
OBJECTS = maind.o gutsf.o gutsp.o misc.o boundary.o part_init.o initial.o

hybrid:	$(OBJECTS) 
	$(F77) -o hybrid $(OBJECTS) 

clean:
	rm *.o hybrid *.out

maind.o:maind.f $(INCLUDE);$(F77) -c maind.f
gutsf.o:gutsf.f $(INCLUDE);$(F77) -c gutsf.f
gutsp.o:gutsp.f $(INCLUDE);$(F77) -c gutsp.f
misc.o:misc.f $(INCLUDE);$(F77) -c misc.f
boundary.o:boundary.f $(INCLUDE);$(F77) -c boundary.f
part_init.o:part_init.f $(INCLUDE);$(F77) -c part_init.f
initial.o:initial.f $(INCLUDE);$(F77) -c initial.f



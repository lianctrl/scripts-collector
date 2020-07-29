FC=gfortran
#FC=g95
#STDF= -std=F
FFLAGS =

TARGET= ../bin/cda 
OBJECTS= cda.o \
         kinds.o \
         paths.o \
         cubes.o

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJECTS)

cda.o: cda.f95 kinds.o cubes.o paths.o
cubes.o: cubes.f95 kinds.o
paths.o: paths.f95 kinds.o
kinds.o: kinds.f95

.SUFFIXES: .f95 .o
.f95.o:
	$(FC) $(STDF) $(FFLAGS) -c -I. -o $*.o $*.f95

clean:
	$(RM) $(TARGET) *.mod *.o

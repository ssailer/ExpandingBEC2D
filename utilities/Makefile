include ../Makefile.inc 

OBJLIBS	= ../libbh3utilities.a
OBJS	= complexgrid.o cudacomplexgrid.o bh3binaryfile.o gauss_random.o realgrid.o cudarealgrid.o

CXXFLAGS	= $(PRJCXXFLAGS) -I. -I../propagators -I../startgrids -I../observables
CUFLAGS		= $(PRJCUFLAGS) -I. -I../propagators -I../startgrids -I../observables



all: $(OBJLIBS)
	echo hello $(CXXFLAGS)

../libbh3utilities.a: $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libbh3utilities.a $?
	$(AR) $(ARFLAGS) rv ../libbh3utilities.a $?
	$(ECHO) $(RANLIB) ../libbh3utilities.a
	$(RANLIB) ../libbh3utilities.a

clean:
	$(ECHO) cleaning up in bh3utilities
	-$(RM) -f $(OBJS)


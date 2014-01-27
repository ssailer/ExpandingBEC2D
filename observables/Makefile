include ../Makefile.inc 

OBJLIBS	= ../libbh3observables.a
OBJS	= bh3observables.o

CXXFLAGS	= $(PRJCXXFLAGS) -I../utilities -I../propagators -I../startgrids -I.
CUFLAGS		= $(PRJCUFLAGS) -I../utilities -I../propagators -I../startgrids -I.

all: $(OBJLIBS)

../libbh3observables.a: $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libbh3observables.a $?
	$(AR) $(ARFLAGS) rv ../libbh3observables.a $?
	$(ECHO) $(RANLIB) ../libbh3observables.a
	$(RANLIB) ../libbh3observables.a

clean:
	$(ECHO) cleaning up in bh3observables
	-$(RM) -f $(OBJS)


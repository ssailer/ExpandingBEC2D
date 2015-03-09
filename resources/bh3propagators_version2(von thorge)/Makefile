include ../Makefile.inc 

OBJLIBS	= ../libbh3propagators.a
OBJS	= bh3propagator.o bh3cudapropagator.o 

CXXFLAGS	= $(PRJCXXFLAGS) -I../bh3utilities -I. -I../bh3startgrids -I../bh3observables
CUFLAGS		= $(PRJCUFLAGS) -I../bh3utilities -I. -I../bh3startgrids -I../bh3observables

all: $(OBJLIBS)

../libbh3propagators.a: $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libbh3propagators.a $?
	$(AR) $(ARFLAGS) rv ../libbh3propagators.a $?
	$(ECHO) $(RANLIB) ../libbh3propagators.a
	$(RANLIB) ../libbh3propagators.a

clean:
	$(ECHO) cleaning up in bh3propagators
	-$(RM) -f $(OBJS)


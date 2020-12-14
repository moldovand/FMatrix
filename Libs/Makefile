CC	= 	g++
CFLAGS	= 	-O3

OBJS	= 	VectorN.o Vector.o \
		MatrixMN.o Matrix.o FileScan.o \
		ImagePosition.o SpacePosition.o \
		alg_eigen.o alg_lu.o alg_matutil.o

TOBJS	=	tensor.o transform.o

libVMSC.a : $(OBJS)
	ar r $@ $(OBJS)
	- ranlib $@

libtool.a : $(TOBJS)
	ar r $@ $(TOBJS)
	- ranlib $@

clean :
	rm -f $(OBJS) $(TOBJS)

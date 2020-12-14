CC	= g++
CFLAGS	= -g -O3 -Wall

SRC	= main.c renormal.c pdp.c correct.c extra.c filedeal.c

fundamental : $(SRC) Libs/libVMSC.a Libs/libtool.a
	$(CC) $(CFLAGS) -ILibs $(SRC) -lm -LLibs -lVMSC -ltool -o $@

Libs/libVMSC.a :
	(cd Libs; make libVMSC.a; make clean)

Libs/libtool.a :
	(cd Libs; make libtool.a; make clean)

clean:
	/bin/rm -rf fundamental*

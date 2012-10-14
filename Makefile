CC = g++
CFLAGS = -O2
OBJS = rmfeos.o
LIBS = -lm -lgsl -lgslcblas

rmfeos : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o rmfeos $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -c $<
	
clean:
	rm -f *.o rmfeos

tar:
	tar -zcvf rmfeos$(shell date +%d%m%y).tar.gz *.c *.h Makefile
	

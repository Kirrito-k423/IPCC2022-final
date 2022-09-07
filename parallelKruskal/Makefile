OPENMP="y"
include common/parallelDefs

BENCH = MST
OBJS = MST.o

# include common/MakeBenchLink
MST.o: MST.C
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c MST.C

clean:
	rm -f *.o

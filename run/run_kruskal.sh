source /opt/intel/oneapi/setvars.sh
/opt/intel/oneapi/mpi/2021.6.0/bin/mpicxx ../feGRASS/kruskal_mpi.cpp -o kruskal
mpirun -np 2 ./kruskal byn1.mtx
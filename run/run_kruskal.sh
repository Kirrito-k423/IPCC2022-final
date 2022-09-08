source /opt/intel/oneapi/setvars.sh
rm -rf kruskal
/opt/intel/oneapi/mpi/2021.6.0/bin/mpicxx -g  ../feGRASS/kruskal_mpi.cpp  -o kruskal
# 因为使用的是二叉树归约，进程数需要为2的指数次方
/opt/intel/oneapi/mpi/2021.6.0/bin/mpirun -np 4 ./kruskal g-3989-158400.mtx 2>&1 | tee log.txt

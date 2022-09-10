#!/bin/bash
###
 # @Descripttion: 
 # @version: 
 # @Author: Shaojie Tan
 # @Date: 2022-08-26 15:39:02
 # @LastEditors: Shaojie Tan
 # @LastEditTime: 2022-09-04 11:48:50
### 
if [ $# -eq 0 ]; then
    case=0
elif [[ $# -eq 1 || $# -eq 2 ]]; then
    case=$1
else
    echo "Usage: $0 [test case] [is Debug]"
    exit 0
fi
echo "Test case ${case}"

set -o xtrace   #开启命令显示
export OMP_PROC_BIND=close;export OMP_PLACES=cores

#check if build exist
if [ ! -d ../build/bin ];then
    mkdir -p ../build/bin
fi
if [ ! -d ../build/obj ];then
    mkdir -p ../build/obj
fi

timestamp=`date +"%Y-%m-%d~%H"`
LOG=run_case${case}.$timestamp.log

#get commit
git log -1 > $LOG

#build
if [ $# -eq 2 ]; then
    if [ $2 == "debug" ]; then
        make debugMpi -C ..
    else
        make timeMpi -C ..
    fi
else
    make mpi -C ..
fi

#check if make success
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi

coord_file_list=(./byn.mtx ./byn1.mtx ./g-3989-158400.mtx ./g-6997-83688.mtx ./ks2010-conv.mtx ./ibmpg5.mtx ./ibmpg6.mtx ./ibmpg3.mtx ./ibmpg4.mtx)
ref_file_list=(refer-byn.txt refer-byn1.txt refer-3989-158400.txt refer-6997-83688.txt refer-ks2010.txt refer-ibmpg5.txt refer-ibmpg6.txt refer-ibmpg3.txt refer-ibmpg4.txt)
#run
stdbuf --output=L mpirun -np 2 ../build/bin/main ${coord_file_list[${case}]} 2>&1 |tee -a $LOG

python3 check.py result.txt ${ref_file_list[${case}]} 2>&1 |tee -a $LOG
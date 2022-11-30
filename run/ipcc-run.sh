#!/bin/bash
if [ $# -eq 0 ]; then
    case=0
elif [[ $# -eq 1 || $# -eq 2 ]]; then
    case=$1
else
    echo "Usage: $0 [test case] [is Debug]"
    exit 0
fi
echo "Test case ${case}"

source /public1/soft/modules/module.sh
module load gcc/10.2.0
set -o xtrace   #开启命令显示
export OMP_PROC_BIND=close;export OMP_PLACES=cores
#export OMP_NUM_THREADS=32

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
        make debugPrint -C ..
    else
        make timePrint -C ..
    fi
else
    make -C ..
fi

#check if make success
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi

coord_file_list=(./g-3989-158400.mtx ./g-6997-83688.mtx ./g-15991-606980.mtx)
ref_file_list=(refer-3989-158400.txt refer-6997-83688.txt refer-15991-606980.txt)
#run
stdbuf --output=L srun -p IPCC -N 1 -c 64 -t 2 ../build/bin/main ${coord_file_list[${case}]} 2>&1 |tee -a $LOG

err_count=`diff result.txt ${ref_file_list[${case}]} | wc -l`
if [ $err_count -eq 0 ]; then
    echo "Test case ${case} passed" |tee -a $LOG
else
    echo "Test case ${case} failed" |tee -a $LOG
fi
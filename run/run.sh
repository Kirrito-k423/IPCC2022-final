#!/bin/bash
###
 # @Descripttion: 
 # @version: 
 # @Author: Shaojie Tan
 # @Date: 2022-08-26 15:39:02
 # @LastEditors: Shaojie Tan
 # @LastEditTime: 2022-08-26 15:43:46
### 
if [ $# -eq 0 ]; then
    case=0
elif [ $# -eq 1 ]; then
    case=$1
else
    echo "Usage: $0 [test case]"
    exit 0
fi
echo "Test case${case}"

set -o xtrace   #开启命令显示
export OMP_PROC_BIND=close;export OMP_PLACES=cores

#check if build exist
if [ ! -d ../build/bin ];then
    mkdir -p ../build/bin
fi

timestamp=$(git log --pretty=format:"%cd" -1 --date=format:'%Y-%m-%d_%H:%M:%S')
LOG=run.$timestamp.log

#get commit
git log -1 > $LOG

#build
make -C ..

#check if make success
if [ $? -ne 0 ]; then
    echo "Make failed"
    exit 1
fi

coord_file_list=(./byn.txt ./byn1.txt)
ref_file_list=(refer-byn.txt refer-byn.txt)
#run
stdbuf --output=L ../build/bin/main ${coord_file_list[${case}]} 2>&1 |tee -a $LOG

#check if run success
if [ $? -ne 0 ]; then
    echo "Run failed"
    exit 1
fi
python3 check.py result.txt ${ref_file_list[${case}]} 2>&1 |tee -a $LOG
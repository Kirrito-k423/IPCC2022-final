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
    echo "case为数字：输入和参考文件从列表中选取"
    echo "case为文件路径：输入为指定文件，参考文件为文件名加上'ref-'"
    echo "      ex. 输入为：thupg1.mtx_layer_2，参考文件为：ref-thupg1.mtx_layer_2"
    exit 0
fi
echo "Test case ${case}"

set -o xtrace   #开启命令显示
export OMP_PROC_BIND=close;export OMP_PLACES=cores

#check if build exist
if [ ! -d ../build/bin ];then
    mkdir -p ../build/bin
fi

timestamp=`date +"%Y-%m-%d~%H"`
if [[ $case =~ ^[0-9]+$ ]]; then
    LOG=run_case${case}.$timestamp.log
else
    basename=$(python <<EOF
var="$case"
var=var.split('/')[-1].split('.')[0]
print(var)
EOF
) 
    LOG=run_case_${basename}.$timestamp.log
fi

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

coord_file_list=(./byn.mtx ./byn1.mtx ./g-3989-158400.mtx ./g-6997-83688.mtx ./ks2010-conv.mtx ./ibmpg5.mtx ./ibmpg6.mtx ./ibmpg3.mtx ./ibmpg4.mtx)
ref_file_list=(refer-byn.txt refer-byn1.txt refer-3989-158400.txt refer-6997-83688.txt refer-ks2010.txt refer-ibmpg5.txt refer-ibmpg6.txt refer-ibmpg3.txt refer-ibmpg4.txt)
#run
if [[ $case =~ ^[0-9]+$ ]]; then
    input=${coord_file_list[${case}]}
    ref=${ref_file_list[${case}]}
else
    input=$case
    ref=$(python <<EOF
var="$input"
var=var.split('/')[-1]
print("ref-"+var)
EOF
)
fi
stdbuf --output=L ../build/bin/main ${input} 2>&1 |tee -a $LOG

python3 check.py result.txt ${ref} 2>&1 |tee -a $LOG
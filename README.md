## 运行方式

在run目录下提供了ipcc-run.sh脚本。该脚本会自动编译代码，运行选择的case，最终通过srun在超算上运行。如下所示：
```bash
srun -p IPCC -N 1 -c 64 -t 2 ../build/bin/main ${coord_file_list[${case}]}
```

1. 将3个case数据和refer参考结果放至run目录下
2. 进入run目录，按照下面命令运行3个case
 - *注：通过`export OMP_NUM_THREADS`设置线程数。case3(g-3989-158400)和case4(g-6997-83688)使用32线程，case5(g-15991-606980)使用64线程*

```bash
export OMP_NUM_THREADS=32 ./ipcc-run.sh 0   # case3
export OMP_NUM_THREADS=32 ./ipcc-run.sh 1   # case4
export OMP_NUM_THREADS=64 ./ipcc-run.sh 2   # case5
```
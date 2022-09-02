# 检查排序后的copy_off_tree_edge
# python check-effect-resist.py case2-copy_off_tree_edge-sort.log edge-copy_off_tree_edge-sort.log
import sys

if len(sys.argv)!=3:
    print("usage: check file file_ref")
    exit(-1)

file = sys.argv[1]
file_ref = sys.argv[2]

def read_pivots_list(filename):
    pivots_list = []
    fp = open(filename, "r")
    for line in fp.readlines(): #读取每一行
        pivots = [int(pivot) for pivot in line.strip().split(' ')[:2]]  #pivot由空格分隔，strip去除行尾回车
        pivots_list.append(pivots)
    fp.close()
    return pivots_list

pivots_list = read_pivots_list(file)
pivots_list_ref = read_pivots_list(file_ref)

#比较
#赛题只说了pviots之间需要有序，没有说明pviots内部，是否需要按索引有序。
for i in range(len(pivots_list)):
    s1 = set(pivots_list[i])    
    s2 = set(pivots_list_ref[i])
    if s1!=s2:
        print(i, s1, s2)
        exit(-1)
print("Check Passed")
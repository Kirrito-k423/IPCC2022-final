# ks2010提供的数据不包含对称重复的边，需要转换
# python ks2010.mtx ks2010-conv.mtx
import sys

if len(sys.argv)!=3:
    print("usage: prog fin fout")
    exit(-1)

file = sys.argv[1]
file_out = sys.argv[2]

def read_edges_swap(filename):
    edge_list = []
    fp = open(filename, "r")
    head = fp.readline().strip().split(' ')
    for line in fp.readlines(): #读取每一行
        split_list = line.strip().split(' ')
        edge_list.append([int(split_list[0]), int(split_list[1])] + split_list[2:])
        edge_list.append([int(split_list[1]), int(split_list[0])] + split_list[2:])
    fp.close()
    return head, edge_list

head_, edge_list = read_edges_swap(file)
head = [int(head_[0]), int(head_[1]), int(head_[2])*2]
edge_list.sort(key = lambda x: x[1])


# for i in range(10):
#     line = edge_list[i]
#     print("%d %d %s"%(line[0], line[1], line[2]))

with open(file_out, "w") as fp:
    fp.write("%d %d %d\n"%(head[0], head[1], head[2]))
    for line in edge_list:
        fp.write("%d %d %s\n"%(line[0], line[1], line[2]))




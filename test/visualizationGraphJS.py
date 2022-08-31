'''
Descripttion: 
version: 
Author: Shaojie Tan
Date: 2022-08-31 19:08:33
LastEditors: Shaojie Tan
LastEditTime: 2022-08-31 19:52:15
'''
import numpy as np
from icecream import ic
# https://observablehq.com/d/10bbfea77cb717b4

f = open("../run/byn1.mtx", "r")
lineList = f.readline() .split()
M = int(lineList[0])
edgeNum = lineList[2]
ic(M, edgeNum) 
# 创建一个 3x4 的数组且所有值全为 0
outMatrix = np.zeros((M, M), dtype=int)
for line in f.readlines():
    lineList = line.split()
    ic(lineList)
    point1=int(lineList[0])-1
    point2=int(lineList[1])-1
    w = int(lineList[2])
    outMatrix[point1][point2]=w
    outMatrix[point2][point1]=w

ic(outMatrix)



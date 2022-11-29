<!--
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-26 15:21:38
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-08-27 16:31:51
-->
# IPCC2022-final

## 运行

### 小例子
```
cd run
./run.sh 1
./run.sh 1 0        # 打印一些基本时间信息
./run.sh 1 debug # 打印debug信息
```

## 优化过程

- v0.1 （baseline）
- v0.2 (简单OMP并行)
    - 对恢复边阶段粗粒度并行
- v0.3（使用LCA计算生成树上等效电阻）
- v0.4（使用邻接表）
    - 计算等效电阻时产生了生成树的邻接表。可以用于恢复边阶段的BFS
- v0.5（使用hash表加速非树边相似性判断）
    - map: (u, v) -> 非树边索引
- v0.6（粗粒度OMP部分实现改用similar_list）
    - 原本粗粒度部分通过数组标记相似边，使用vector避免分配巨大的数组。
- v0.7（细粒度、粗粒度两阶段OMP）
- v1.0 （使用并查集优化串行kruskal）
    - 原本kruskal的实现，也使用了并查集，但是每次unite操作时，会进行路径压缩，将所有点parent设置为另一方根节点。
- v1.1（每个顶点维护一个相邻、非树边map）
    - G_adja(u):  map: (v) -> 非树边索引
    - 比使用一个大map快许多，值得考虑什么时候单独map更快？
- v1.2（计算等效电阻部分使用omp并行）
- v1.3（使用结构体代替vector<double>）
- v1.4（替换并行排序）
    - 自己实现的归并排序
- v1.5（相似性判断交换循环顺序）
- v1.6  (优化掉delta计算mark的使用)
- v1.7 （树上并行BFS)
    - 对thupg大例子效果明显（大例子BFS是热点）
- v2.1（PSRS排序）
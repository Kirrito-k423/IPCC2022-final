/*
 * @Descripttion:
 * @version:
 * @Author: Shaojie Tan
 * @Date: 2022-08-27 15:58:17
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-05 17:19:46
 */
#ifndef _GLOBAL_H
#define _GLOBAL_H
// #include <omp.h>
// #include "Eigen/Dense"
// #include "Eigen/LU"
#include <math.h>
#include <string.h>
#include <sys/time.h>

// #include <algorithm>
#include <parallel/algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include <unordered_set>

extern int NUM_THREADS;         //全局线程数，在main.cpp中定义。
                                //通过export OMP_NUM_THREADS设置，通过omp_get_max_threads()获取

#define CREATE_ADJA_THREADS     NUM_THREADS     //等效电阻部分构造生成树的邻接表
#define SORT1_THREADS           NUM_THREADS     //kruskal算法，对图G的边按有效权重排序
#define SORT2_THREADS           NUM_THREADS     //恢复边阶段之前，对非树边对等效电阻降序排序
#define SORT3_THREADS           NUM_THREADS     //细粒度排除相似边中，对bfs_process排序
#define TREE_BFS_THREADS        NUM_THREADS     //树上并行进行beta_BFS的线程数（对大例子有效，小例子无效）
#define TREE_BFS_FACTOR         3               //在生成树上进行BFS，并行第TREE_BFS_LAYER层
#define first_step_OMP_percentage 0.00          //细粒度并行部分，解决边数的占比（对大例子有效，小例子无效）
// #define task_pool_size 512                   //变成由M L 确定的全局变量
#define avail_percent 0.92
#define search_block_size_start 512
#define OFFSET 10

#define cut_similarity_range 17

#define next_range 128 // DEBUG_PRINT

// using namespace Eigen;
using namespace std;

#ifdef TIME
#define TIME_PRINT(fmt, args...) fprintf(stderr, fmt, ##args)
#define printTime(s)                                                                                                                              \
    {                                                                                                                                             \
        gettimeofday(&endTime, NULL);                                                                                                             \
        fprintf(stderr, "%-40s: took %.6f ms\n", s, (endTime.tv_sec - startTime.tv_sec) * 1000 + (endTime.tv_usec - startTime.tv_usec) / 1000.0); \
        gettimeofday(&startTime, NULL);                                                                                                           \
    }
#else
#define TIME_PRINT(fmt, args...) /* Don't do anything in release builds */
#define printTime(s)
#endif

#ifdef DEBUG
#define PIntWithName(var) printf("%s\t%d\n", #var, var)
#define DEBUG_PRINT(fmt, args...) fprintf(stderr, fmt, ##args)
#define OMP_TIME_PRINT(fmt, args...) fprintf(stderr, fmt, ##args)
#else
#define PIntWithName(var) /* Don't do anything in release builds */
#define DEBUG_PRINT(fmt, args...)
#define OMP_TIME_PRINT(fmt, args...)
#endif

struct edge{
    int u, v;
    double eff_w, w;
    edge(){}
    edge(int u, int v, double eff_w, double w) :u(u), v(v), eff_w(eff_w), w(w){}
    bool operator<(const edge& e) const noexcept {
        if(eff_w == e.eff_w){
            if(v == e.v){
                return u > e.u;
            }
            return v > e.v;
        }
        return eff_w < e.eff_w;
    }
    bool operator==(const edge& e) const noexcept {
        return u==e.u && v==e.v && eff_w==e.eff_w && w==e.w;
    }
};
typedef struct edge edge_t;
int cmp_by_effw(const void *a, const void *b);
struct adj_node{
    int u;
    double w;
};
typedef struct adj_node node_t;

// global value
extern int M;
extern int N;
extern int L;
extern int largest_volume_point;

extern vector<vector<int>> adja_list;
extern double *dis;
extern int *parent;
extern int *no_weight_dis;
extern double fg_similarity_time[2];                   // 伪逆， 循环总时间， 循环内三部分时间

// recover_off_edges.cpp
int calculate_beta(int i, int j);
void beta_BFS(int beta, std::vector<int> &queue, int root);
void beta_BFS_p(int beta, std::vector<int> &queue, int root);
void adjust_similarity_tree(std::vector<int> &bfs_process1, std::vector<int> &bfs_process2,
                            vector<int> &similar_list, vector<vector<std::pair<int, int>>> &G_adja);

void fg_adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2,
                               int *similarity_tree, vector<vector<std::pair<int, int>>> &G_adja);
void check_next_range_similarity_tree(int i, int *similarity_tree, int total_range);
void merge_thread_similarity_tree(int i, int similarity_tree_length, int *similarity_tree, int *thread_similarity_tree_address);

// effect_resistance.cpp
void caculate_resistance(vector<edge_t> &spanning_tree, vector<edge_t> &off_tree_edge, vector<edge_t> &copy_off_tree_edge);
void write_edge(vector<edge_t> &edge, const char *file);
void write_bfs_process(vector<int> &bfs_process1, const char *file);
int get_LCA(int i, int j, int *parent, int *no_weight_dis);

void print_M1_Array(string name, int *toPrint);
void printStack(string name, stack<int> toPrint);
int get_task_pool_size(int total_num);

void kruscal(vector<edge_t> &edge_matrix, vector<edge_t> &spanning_tree);
bool compare(const edge_t &a, const edge_t &b);
int cmp(const void *a, const void *b);
#include "p_mergesort.hpp"
void psrs(vector<edge_t> &arr, int p, comparison_fn_t cmp);
#endif

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
#include <stack>
#include <map>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "common/graph.h"
#include "parlay/sequence.h"
#include "MST.h"



#define NUM_THREADS 32
#define first_step_OMP_percentage 0.02  //第一部分OMP的解决边数的占比 case2 3 0.01更快
// #define task_pool_size 512           //变成由M L 确定的全局变量
#define avail_percent 0.92
#define search_block_size_start 512
#define offset 10

#define cut_similarity_range 3


#define next_range 128 //DEBUG_PRINT

// using namespace Eigen;
using namespace std;

#ifdef TIME
#define TIME_PRINT(fmt, args...)    fprintf(stderr, fmt, ## args)
#else
#define TIME_PRINT(fmt, args...)    /* Don't do anything in release builds */
#endif

#ifdef DEBUG
#define PIntWithName(var)  printf("%s\t%d\n", #var, var) 
#define DEBUG_PRINT(fmt, args...)    fprintf(stderr, fmt, ## args)
#define OMP_TIME_PRINT(fmt, args...)    fprintf(stderr, fmt, ## args)
#else
#define PIntWithName(var)    /* Don't do anything in release builds */
#define DEBUG_PRINT(fmt, args...)  
#define OMP_TIME_PRINT(fmt, args...)
#endif

#define printTime(j) {\
    gettimeofday(&endTime, NULL);\
    TIME_PRINT(j, (endTime.tv_sec-startTime.tv_sec)*1000+(endTime.tv_usec-startTime.tv_usec)/1000.0);\
    gettimeofday(&startTime, NULL);\
}

//global value
extern int M;
extern int N;
extern int L;
extern int largest_volume_point;

extern vector<vector<vector<double>>> adja_list;
extern double *dis;
extern int *parent;
extern int *no_weight_dis;

//recover_off_edges.cpp
int calculate_beta(int i, int j);
void beta_BFS(int beta, std::vector<int> &queue, int root);
void adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, vector<vector<double>> &copy_off_tree_edge);
void adjust_similarity_tree(std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            vector<int> &similar_list, map<uint64_t, uint32_t> &copy_off_tree_edge);
void fg_adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, map<uint64_t, uint32_t> &off_tree_edge_map);
void check_next_range_similarity_tree(int i, int *similarity_tree, int total_range);
void merge_thread_similarity_tree(int i, int similarity_tree_length, int * similarity_tree, int *thread_similarity_tree_address);

//effect_resistance.cpp
void caculate_resistance(vector<vector<double>> &spanning_tree, vector<vector<double>> &off_tree_edge, vector<vector<double>> &copy_off_tree_edge);
void write_edge(vector<vector<double>> &edge, const char *file);
int get_LCA(int i, int j, int *parent, int *no_weight_dis);

void print_M1_Array(string name,int * toPrint);
void printStack(string name, stack<int> toPrint);
int get_task_pool_size(int total_num);

template <class intV, class Weight>
wghEdgeArray<intV,Weight> read_WghEdgeArray_From_edge_matrix(vector<vector<double>> &W) {
  using WE = wghEdge<intV,Weight>;
  // parlay::sequence<char> S = readStringFromFile(fname);
  // parlay::sequence<char*> W = stringToWords(S);
  // if (W[0] != WghEdgeArrayHeader) {
  //   cout << "Bad input file" << endl;
  //   abort();
  // }
  long n = W.size();
  cout << W.size() << "\n" 
  << n << "\n" 
  << W[0][0] << "\n"
  << W[0][1] << "\n" 
  << W[0][2] << "\n" 
  << W[0][3] << "\n" 
//   << W[0][4] << "\n" 
//   << W[0][5] << "\n"
  << W[1][0] << "\n"
  << W[1][1] << "\n" 
  << W[1][2] << "\n" 
  << W[1][3] << "\n" 
  << std::endl;
  auto E = parlay::tabulate(n, [&] (size_t i) -> WE {
return WE((long int)(W[i][0]),
    (long int)(W[i][1]),
    (Weight) (-W[i][2]));});

  auto mon = parlay::make_monoid([&] (WE a, WE b) {
return WE(std::max(a.u, b.u), std::max(a.v, b.v), 0);},
    WE(0,0,0));
  auto r = parlay::reduce(E, mon);

  return wghEdgeArray<intV,Weight>(std::move(E), max<intV>(r.u, r.v) + 1);
}
#endif
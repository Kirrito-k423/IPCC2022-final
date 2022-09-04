/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-27 15:58:17
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-03 09:37:15
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

#define NUM_THREADS 32
#define task_pool_size 512
// enum task_divide_mode {sequential, distracted}
// #define mode sequential
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
                            vector<int> &similar_list, map<uint32_t, uint16_t> &copy_off_tree_edge);
void check_next_range_similarity_tree(int i, int *similarity_tree, int total_range);
void merge_thread_similarity_tree(int i, int similarity_tree_length, int * similarity_tree, int *thread_similarity_tree_address);

//effect_resistance.cpp
void caculate_resistance(vector<vector<double>> &spanning_tree, vector<vector<double>> &off_tree_edge, vector<vector<double>> &copy_off_tree_edge);
void write_edge(vector<vector<double>> &edge, const char *file);
int get_LCA(int i, int j, int *parent, int *no_weight_dis);

void print_M1_Array(string name,int * toPrint);
void printStack(string name, stack<int> toPrint);
#endif
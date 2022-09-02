/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-27 15:58:17
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-02 12:46:11
 */
#ifndef _GLOBAL_H
#define _GLOBAL_H
// #include <omp.h>
#include "Eigen/Dense"
#include "Eigen/LU"
#include <stack>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <sys/time.h>

#define NUM_THREADS 16
#define task_pool_size 16
// enum task_divide_mode {sequential, distracted}
// #define mode sequential
#define cut_similarity_range 3
#define next_range 128 //DEBUG_PRINT

using namespace Eigen;
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


//global value
extern int M;
extern int N;
extern int L;


int calculate_belta(int i, MatrixXd &LG, int largest_volume_point, int edge_point1, int edge_point2);
void belta_BFS(int belta, MatrixXd &LG, std::vector<int> &candidate_point_set, int edge_point);
void adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, vector<vector<double>> &copy_off_tree_edge);
void check_next_range_similarity_tree(int i, int *similarity_tree, int total_range);

void merge_thread_similarity_tree(int i, int similarity_tree_length, int * similarity_tree, int *thread_similarity_tree_address);
void test_LCA_find_update(MatrixXd &LG, int largest_volume_point);
void LCA_find(int *find, MatrixXd &LG, int largest_volume_point);
int calculate_belta_from_find(int i, int *find, int edge_point1, int edge_point2);

void print_M1_Array(string name,int * toPrint);
void printStack(string name, stack<int> toPrint);
#endif
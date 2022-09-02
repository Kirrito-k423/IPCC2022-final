/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-27 15:58:17
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-08-30 12:59:24
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


int calculate_belta(int i, MatrixXd *LG, int largest_volume_point, int edge_point1, int edge_point2);
void belta_BFS(int belta, MatrixXd *LG, std::vector<int> *candidate_point_set, int edge_point);
void adjust_similarity_tree(int i, std::vector<int> *bfs_process1, std::vector<int> *bfs_process2 ,\
                            int *similarity_tree, vector<vector<double>> *copy_off_tree_edge);
void check_next_range_similarity_tree(int i, int *similarity_tree, int total_range);

void caculate_resistance(vector<vector<double>> &spanning_tree, vector<vector<double>> &off_tree_edge, vector<vector<double>> &copy_off_tree_edge);
void write_edge(vector<vector<double>> &edge, const char *file);
// void write_edge(vector<vector<double>> &edge, string file);
#endif
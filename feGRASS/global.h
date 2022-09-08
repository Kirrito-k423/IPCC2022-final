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
#include "mpi.h"
#include<signal.h>

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
void kruscal(vector<vector<double>> &edge_matrix, vector<vector<double>> &spanning_tree);
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

// parallel kruskal using mpi
typedef enum { FALSE,
               TRUE } boolean;

void parse_input(vector<vector<double>> &edges, int rank, int number_of_edges, int vertices_per_process, int edges_per_process);
void parse_edge_list_input();
int compare_edges(const void *a, const void *b);
boolean should_send();
int get_merge_partner_rank(boolean should_send, int rank, int merge_iteration);
void merge(int rank, int number_of_vertices, vector<vector<double>> &local_msf_edges, vector<vector<double>> &recv_msf_edges, vector<vector<double>> &merged_msf_edges);
void send_recieve_local_msf(int rank, int count, vector<vector<double>> &local_msf_edges, vector<vector<double>> & recv_msf_edges);
void merge_msf(int num_of_processors, int rank, int &merge_iteration);
void sig_handler(int sig) {
	std::cerr << "Slave signal received " << sig << std::endl;
	while (1);
}
typedef struct node_t { /* Union-Find data structure */
    struct node_t *parent;
    int depth;
} u_node;

u_node *uf_set; /* Array indicating which set does vertex belong to
                   (used in Union-Find algorithm) */

void uf_make(int number_of_vertices);
u_node *uf_find(u_node *a);
void uf_union(u_node *a, u_node *b);
#endif
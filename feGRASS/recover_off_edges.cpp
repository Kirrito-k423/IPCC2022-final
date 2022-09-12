/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-29 19:59:51
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-05 16:19:17
 */
#include "global.h"

int calculate_beta(int i, int j){
    int LCA_point = get_LCA(i-1, j-1, parent, no_weight_dis);
    int d1 = no_weight_dis[i-1] - no_weight_dis[LCA_point];
    int d2 = no_weight_dis[j-1] - no_weight_dis[LCA_point];
    return d1 < d2 ? d1 : d2;
}

void beta_BFS(int beta, std::vector<int> &queue, int root){
    queue.push_back(root);
    //use zero to cut the near layer
    queue.push_back(0);
    int * mark = (int * )malloc(M * sizeof(int));
    memset(mark, 0, M *sizeof(int));
    mark[root-1]=1;
    int layer=0;
    for (int j=0;j<M;j++) {
        if (layer==beta){
            break;
        }
        if (queue[j]==0){
            queue.push_back(0);
            layer++;
        }
        else{
            int point = queue[j]-1;
            for (int i=0; i<adja_list[point].size(); i++) {
                int search_point = adja_list[point][i].u;
                if(mark[search_point]==0){
                    queue.push_back(search_point+1);
                    mark[search_point] = 1;
                }
            }
        }
    }
    free(mark);
}

// fine_grained 细粒度
void fg_adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            vector<vector<int>> &similar_list, vector<map<int, int>> &G_adja){
    //mark the edge that is similar to the edge which wants to be added
    int MPI_start = mpi_rank * bfs_process1.size() / comm_size;
    int MPI_end = (mpi_rank+1) * bfs_process1.size() / comm_size;
    MPI_DEBUG_PRINT("mpi_rank %d\t start end %d\t%d\n",mpi_rank,MPI_start,MPI_end);

    //dynamic 会产生 大约60000* 60000 次omp 线程创建开销
    #pragma omp parallel for num_threads(NUM_THREADS) collapse(2)
    for (int j=MPI_start; j<MPI_end; j++) {
        for (int k=0; k<bfs_process2.size(); k++) {
            int tid = omp_get_thread_num();
            if (bfs_process2[k]==0 ||bfs_process1[j]==0) {
                continue;
            }
            if (bfs_process1[j]==bfs_process2[k]) {
                continue;
            }
            int u = bfs_process1[j]-1;
            int v = bfs_process2[k]-1;
            if(G_adja[u].count(v)==1){
               similar_list[tid].push_back(G_adja[u].find(v)->second);
            }
        }
    }
}

void adjust_similarity_tree(std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                         vector<int> &similar_list, vector<map<int, int>> &G_adja){
    for (int j=0; j<bfs_process1.size(); j++) {
        if (bfs_process1[j]==0) {
            continue;
        }
        for (int k=0; k<bfs_process2.size(); k++) {
            if (bfs_process2[k]==0) {
                continue;
            }
            if (bfs_process1[j]==bfs_process2[k]) {
                continue;
            }
            int u = bfs_process1[j]-1;
            int v = bfs_process2[k]-1;
            if(G_adja[u].count(v)==1){
               similar_list.push_back(G_adja[u].find(v)->second);
            }
        }
    }
}

void check_next_range_similarity_tree(int i, int *similarity_tree, int total_range){
    int eqaul_zero_num=0;
    for (int j=i+1; j<= i+total_range; j++){
        if(similarity_tree[j]==0){
            eqaul_zero_num++;
        }
    }
    DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t check_next_range\t %d/%d \t%.2f%%\n"\
                ,i,     eqaul_zero_num,     total_range,     100*(double)eqaul_zero_num/total_range);
}

int get_task_pool_size(int total_num){
    int i = search_block_size_start;
    for(; i < total_num; i++){
        int avail_block_num = i * avail_percent;
        int loop_num = total_num / avail_block_num;
        if(i%2==0 && (loop_num+1) * avail_block_num - total_num < OFFSET){ //粗粒度MPI的原因，要是偶数
            break;
        }
    }
    return i;
}
/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-29 19:59:51
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-05 16:19:17
 */
#include "global.h"

void fg_time_print(struct timeval &start, struct timeval &end, int index){
    gettimeofday(&end, NULL);
    double tmp_past_time = (end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000.0;
    fg_similarity_time[index] += tmp_past_time;
    DEBUG_PRINT("copy_off_tree_edge omp_detail %d \t took %f ms\n", index, tmp_past_time);
}

int cmp2(const void *a, const void *b) {
    return *((int *)a) > *((int *)b);
}

int calculate_beta(int i, int j){
    int LCA_point = get_LCA(i-1, j-1, parent, no_weight_dis);
    int d1 = no_weight_dis[i-1] - no_weight_dis[LCA_point];
    int d2 = no_weight_dis[j-1] - no_weight_dis[LCA_point];
    return d1 < d2 ? d1 : d2;
}

void beta_BFS_p(int beta, std::vector<int> &queue, int root){
    /**
     * 顶点索引从0还是1开始：
     * - adja_list从0开始
     * - queue从1开始
    */
    std::vector<int> parent_queue;
    queue.reserve(beta*100);
    parent_queue.reserve(beta*100);

    queue.push_back(root);
    parent_queue.push_back(0);

    if(beta==0) return;

    //对root节点先处理，避免需要判断root的父节点
    vector<int> &adja = adja_list[root-1];
    for(int search_point: adja){
        queue.push_back(search_point+1);
        parent_queue.push_back(root);
    }

    int layer = 1;
    int begin = 1;
    int last = queue.size();

    while(layer != beta){
        for(int i = begin; i < last; i++){
            int point = queue[i];
            int parent = parent_queue[i];   //parent of point
            vector<int> &adja = adja_list[point-1];
            for(int search_point: adja){
                if(parent != search_point+1){
                    queue.push_back(search_point+1);
                    parent_queue.push_back(point);
                }
            }
        }
        begin = last;
        last = queue.size();
        layer++;
    }
}

void beta_BFS(int beta, std::vector<int> &queue, int root){
    /**
     * 顶点索引从0还是1开始：
     * - adja_list从0开始
     * - queue从1开始
    */
    std::vector<int> parent_queue;
    queue.reserve(beta*100);
    parent_queue.reserve(beta*100);

    queue.push_back(root);
    parent_queue.push_back(0);

    if(beta==0) return;

    //对root节点先处理，避免需要判断root的父节点
    vector<int> &adja = adja_list[root-1];
    for(int search_point: adja){
        queue.push_back(search_point+1);
        parent_queue.push_back(root);
    }

    int layer = 1;
    int begin = 1;
    int last = queue.size();

    while(layer != beta){
        for(int i = begin; i < last; i++){
            int point = queue[i];
            int parent = parent_queue[i];   //parent of point
            vector<int> &adja = adja_list[point-1];
            for(int search_point: adja){
                if(parent != search_point+1){
                    queue.push_back(search_point+1);
                    parent_queue.push_back(point);
                }
            }
        }
        begin = last;
        last = queue.size();
        layer++;
    }
}

// fine_grained 细粒度
void fg_adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, vector<vector<std::pair<int, int>>> &G_adja){
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);
    p_mergesort<int>(bfs_process2, 32, cmp2);
    fg_time_print(start_time,end_time,0);
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int j=0; j<bfs_process1.size(); j++) {
        int u = bfs_process1[j]-1;
        for(auto it=G_adja[u].begin(); it!=G_adja[u].end();it++){
            int v = it->first;
            if(std::binary_search(bfs_process2.begin(), bfs_process2.end(), v+1)){  //v+1
               similarity_tree[it->second] = 1;
            }
        }
    }
    fg_time_print(end_time,start_time,1);
}

void adjust_similarity_tree(std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                         vector<int> &similar_list, vector<vector<std::pair<int, int>>> &G_adja){
    quick_sort<int>(bfs_process2, cmp2);
    for (int j=0; j<bfs_process1.size(); j++) {
        int u = bfs_process1[j]-1;
        for(auto it=G_adja[u].begin(); it!=G_adja[u].end();it++){
            int v = it->first;
            if(std::binary_search(bfs_process2.begin(), bfs_process2.end(), v+1)){  //v+1
               similar_list.push_back(it->second);
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
        if((loop_num+1) * avail_block_num - total_num < OFFSET){
            break;
        }
    }
    return i;
}
/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-29 19:59:51
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-05 16:19:17
 */
#include "global.h"

int cmp2(const void *a, const void *b) {
    return *((int *)a) > *((int *)b);
}

int calculate_beta(int i, int j){
    int LCA_point = get_LCA(i-1, j-1, parent, no_weight_dis);
    int d1 = no_weight_dis[i-1] - no_weight_dis[LCA_point];
    int d2 = no_weight_dis[j-1] - no_weight_dis[LCA_point];
    return d1 < d2 ? d1 : d2;
}

void beta_BFS(int beta, std::vector<int> &queue, int root){
    queue.push_back(root);
    int * mark = (int * )malloc(M * sizeof(int));
    memset(mark, 0, M *sizeof(int));
    mark[root-1]=1;
    int layer=0;
    int begin = 0;
    int last = 1;
    while(layer != beta){
        for(int i = begin; i < last; i++){
            int point = queue[i]-1;
            for(int j = 0; j < adja_list[point].size(); j++){
                int search_point = adja_list[point][j].u;
                if(mark[search_point]==0){
                    queue.push_back(search_point+1);
                    mark[search_point] = 1;
                }
            }
        }
        begin = last;
        last = queue.size();
        layer++;
    }
    // printf("for root: %d, ", root);
    // for(int i = 0; i < queue.size(); i++){
    //     if (queue[i] != 0){
    //         printf(" %d ", queue[i]);
    //     }
    // }
    // printf("\n");
    free(mark);
}

// fine_grained 细粒度
void fg_adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, vector<vector<std::pair<int, int>>> &G_adja){
    phmap::flat_hash_set<int> set;
    for(int i=0; i<bfs_process2.size(); i++){
        set.insert(bfs_process2[i]);
    }
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int j=0; j<bfs_process1.size(); j++) {
        int u = bfs_process1[j]-1;
        for(auto it=G_adja[u].begin(); it!=G_adja[u].end();it++){
            int v = it->first;
            if(set.contains(v+1)){
               similarity_tree[it->second] = 1;
            }
        }
    }
}

void adjust_similarity_tree(std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                         vector<int> &similar_list, vector<vector<std::pair<int, int>>> &G_adja){
    // quick_sort<int>(bfs_process2, cmp2);
    phmap::flat_hash_set<int> set;
    for(int i=0; i<bfs_process2.size(); i++){
        set.insert(bfs_process2[i]);
    }
    for (int j=0; j<bfs_process1.size(); j++) {
        int u = bfs_process1[j]-1;
        for(auto it=G_adja[u].begin(); it!=G_adja[u].end();it++){
            int v = it->first;
            if(set.contains(v+1)){
            // if(std::binary_search(bfs_process2.begin(), bfs_process2.end(), v+1)){  //v+1
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
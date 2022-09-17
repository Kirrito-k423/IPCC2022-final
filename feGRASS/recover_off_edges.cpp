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
void fg_adjust_similarity_tree(int i, std::vector<int> &bfs_process1_, std::vector<int> &bfs_process2_ ,\
                            int *similarity_tree, vector<map<int, int>> &G_adja){
    //mark the edge that is similar to the edge which wants to be added
    std::vector<int> bfs_process1;
    bfs_process1.reserve(bfs_process1_.size());
    std::vector<int> bfs_process2;
    bfs_process2.reserve(bfs_process2_.size());
    for (int j=0; j<bfs_process1_.size(); j++) {
        if(bfs_process1_[j]!=0)
            bfs_process1.push_back(bfs_process1_[j]);
    }
    for (int j=0; j<bfs_process2_.size(); j++) {
        if(bfs_process2_[j]!=0)
            bfs_process2.push_back(bfs_process2_[j]);
    }
    long hit, total, p_pos, p_neg;
    hit = p_pos = p_neg = 0;
    total = bfs_process1.size()*bfs_process2.size();
    #pragma omp parallel for num_threads(NUM_THREADS) collapse(2) reduction(+: hit, p_pos, p_neg)
    for (int j=0; j<bfs_process1.size(); j++) {
        for (int k=0; k<bfs_process2.size(); k++) {
            int u = bfs_process1[j]-1;
            int v = bfs_process2[k]-1;

            if (bloom_contain(u, v)) {
                p_pos++;
                if (G_adja[u].count(v) == 1) {
                    hit++;
                    similarity_tree[G_adja[u].find(v)->second] = 1;
                }
            }
            p_neg++;
        }
    }
    DEBUG_PRINT("hit: %ld, total: %ld, rate: %f\n", hit, total, (double)hit/total);
    DEBUG_PRINT("p_pos: %ld p_neg: %ld p_pos_rate: %f\n", p_pos, p_neg, (double)p_pos/total);
}

void adjust_similarity_tree(std::vector<int> &bfs_process1_, std::vector<int> &bfs_process2_ ,\
                         vector<int> &similar_list, vector<map<int, int>> &G_adja){
    std::vector<int> bfs_process1;
    bfs_process1.reserve(bfs_process1_.size());
    std::vector<int> bfs_process2;
    bfs_process2.reserve(bfs_process2_.size());
    for (int j=0; j<bfs_process1_.size(); j++) {
        if(bfs_process1_[j]!=0)
            bfs_process1.push_back(bfs_process1_[j]);
    }
    for (int j=0; j<bfs_process2_.size(); j++) {
        if(bfs_process2_[j]!=0)
            bfs_process2.push_back(bfs_process2_[j]);
    }
    for (int j=0; j<bfs_process1.size(); j++) {
        for (int k=0; k<bfs_process2.size(); k++) {
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
        if((loop_num+1) * avail_block_num - total_num < OFFSET){
            break;
        }
    }
    return i;
}
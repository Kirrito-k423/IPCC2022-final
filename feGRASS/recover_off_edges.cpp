/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-29 19:59:51
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-09-04 19:14:31
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
    int mark[M];
    memset(mark, 0, sizeof(mark));
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
                int search_point = adja_list[point][i][0];
                if(mark[search_point]==0){
                    queue.push_back(search_point+1);
                    mark[search_point] = 1;
                }
            }
        }
    }
}

void adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, vector<vector<double>> &copy_off_tree_edge){
    //mark the edge that is similar to the edge which wants to be added
    int point_pair=0;
    int hit_num=0;
    int avail_hit=0;

    int hit_cut_num=0;
    int avail_cut_hit=0;

    int hit_next_num=0;
    int avail_next_hit=0;
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) collapse(2)
    for (int j=0; j<bfs_process1.size(); j++) {
        for (int k=0; k<bfs_process2.size(); k++) {
            if (bfs_process2[k]==0 ||bfs_process1[j]==0) {
                continue;
            }
            if (bfs_process1[j]==bfs_process2[k]) {
                continue;
            }
            point_pair++;
            for (int z=i; z<copy_off_tree_edge.size()/cut_similarity_range; z++) { // 余下的off_edge里，如果该边的两点，有一点在两个bfs的点集里，则该边视作similar
                if ((copy_off_tree_edge[z][0]==bfs_process1[j]&&
                    copy_off_tree_edge[z][1]==bfs_process2[k]) ||
                     (copy_off_tree_edge[z][0]==bfs_process2[k]&&
                        copy_off_tree_edge[z][1]==bfs_process1[j])){
                    hit_num++;
                    if(similarity_tree[z]==0)
                        avail_hit++;
                    if(z<i+next_range){
                        hit_next_num++;
                        if(similarity_tree[z]==0)
                            avail_next_hit++;
                    }
                    if(z<copy_off_tree_edge.size()/cut_similarity_range){
                        hit_cut_num++;
                        if(similarity_tree[z]==0)
                            avail_cut_hit++;
                    }
                    similarity_tree[z]=1;
                }
            }
        }
    }
    // DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t hit\t %d \t avail\t %d\t total\t %d\n",i,hit_num,avail_hit,point_pair);
    // DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t hit\t %.2f%% \t avail\t %.2f%%\n",i,100*(double)hit_num/point_pair,100*(double)avail_hit/point_pair);
    // DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t hit\t %d \t avail\t %d \tnext\n",i,hit_next_num,avail_next_hit);
    // DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t hit\t %.2f%% \t avail\t %.2f%%\n",i,100*(double)hit_next_num/point_pair,100*(double)avail_next_hit/point_pair);
    // DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t hit\t %d \t avail\t %d \tcut\n",i,hit_cut_num,avail_cut_hit);
    // DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t hit\t %.2f%% \t avail\t %.2f%%\n",i,100*(double)hit_cut_num/point_pair,100*(double)avail_cut_hit/point_pair);

}


void adjust_similarity_tree(int i, std::vector<int> &bfs_process1, std::vector<int> &bfs_process2 ,\
                            int *similarity_tree, map<uint32_t, uint16_t> &off_tree_edge_map){
    //mark the edge that is similar to the edge which wants to be added
    // int hit_num=0;
    // int avail_hit=0;

    // int hit_cut_num=0;
    // int avail_cut_hit=0;

    // int hit_next_num=0;
    // int avail_next_hit=0;

    //unordered_map 会慢一倍
    map<uint32_t, uint16_t>::iterator endIter = off_tree_edge_map.end();
    //dynamic 会产生 大约60000* 60000 次omp 线程创建开销
    #pragma omp parallel for num_threads(NUM_THREADS) collapse(2)
    for (int j=0; j<bfs_process1.size(); j++) {
        for (int k=0; k<bfs_process2.size(); k++) {
            if (bfs_process2[k]==0 ||bfs_process1[j]==0) {
                continue;
            }
            if (bfs_process1[j]==bfs_process2[k]) {
                continue;
            }
            uint32_t key = (uint32_t(bfs_process1[j]) << 16) | uint32_t(bfs_process2[k]);
            // DEBUG_PRINT("key1 = %x, key2 = %x\n", key1, key2);
            //map<uint32_t, uint16_t>::iterator it;
            map<uint32_t, uint16_t>::iterator iter = off_tree_edge_map.find(key);
            if (iter != endIter) {
                // DEBUG_PRINT("edge index: %d\n", uint16_t(off_tree_edge_map.find(key)->second));
                similarity_tree[uint16_t(iter->second)] = 1;
                // DEBUG_PRINT("hash_edges: node %x, %x, %d\n", uint32_t(bfs_process1[j]), uint32_t(bfs_process2[k]), off_tree_edge_map.find(key)->second);
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

void merge_thread_similarity_tree(int i, int similarity_tree_length, int * similarity_tree, int *thread_similarity_tree_address){
    for(; i < similarity_tree_length; i++){
        if(similarity_tree[i]==0 && thread_similarity_tree_address[i]==1){
            similarity_tree[i]=1;
        }
    }
}
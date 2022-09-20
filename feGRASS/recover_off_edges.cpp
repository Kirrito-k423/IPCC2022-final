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
    queue.push_back(root);
    if(beta==0) return;

    std::vector<int> pqueue;    //记录queue中点对应的parent
    queue.reserve(beta*100);
    pqueue.reserve(beta*100);
    pqueue.push_back(0);

    const int pre_layer = 10;   //先遍历若干层，直到一层有足够多点时，再omp并行

    int layer = 0;
    int begin = 0;
    int last = 1;
    while(layer != beta && layer < pre_layer){
        for(int i = begin; i < last; i++){
            int point = queue[i];
            int parent = pqueue[i];   //parent of point
            for(int search_point: adja_list[point-1]){
                if(parent != search_point+1){
                    queue.push_back(search_point+1);
                    pqueue.push_back(point);
                }
            }
        }
        begin = last;
        last = queue.size();
        layer++;
    }

    if(layer==beta) return;     //如果beta <= pre_layer，则不必再并行了

    const int p = 32;
    vector<vector<int>> queue_(p);      //每个线程各自的队列
    vector<vector<int>> pqueue_(p);

    //任务划分，线程间点的个数相差小于1
    const int n = last - begin;

    int blk_size[p];
    int offset[p + 1]; memset(offset, 0, sizeof(offset)); offset[p] = n;
    for (int i = 0; i < p; i++) {
        blk_size[i] = n / p;
        if (i < n % p) { // r*(q+1) + (p-r)*q
            blk_size[i] += 1;
        }
    }

    #pragma omp parallel num_threads(p)
    {
        int tid = omp_get_thread_num();

        for (int i = 0; i < tid; i++) {
            offset[tid] += blk_size[i];
        }

        queue_[tid].resize(blk_size[tid]);
        pqueue_[tid].resize(blk_size[tid]);
        mempcpy(queue_[tid].data(), queue.data() + offset[tid], blk_size[tid] * sizeof(int));
        mempcpy(pqueue_[tid].data(), pqueue.data() + offset[tid], blk_size[tid] * sizeof(int));

        //各个线程在各自子树上BFS
        int layer_ = layer;
        int begin_ = 0;
        int last_ = blk_size[tid];
        while(layer_ != beta){
            for(int i = begin_; i < last_; i++){
                int point = queue_[tid][i];
                int parent = pqueue_[tid][i];   //parent of point
                for(int search_point: adja_list[point-1]){
                    if(parent != search_point+1){
                        queue_[tid].push_back(search_point+1);
                        pqueue_[tid].push_back(point);
                    }
                }
            }
            begin_ = last_;
            last_ = queue_[tid].size();
            layer_++;
        }
    }

    //合并
    int cnt = 0;
    for(int i=0; i<p; i++){
        // memcpy(queue.data() + cnt, queue_[i].data(), queue_[i].size()*sizeof(int));
        // cnt += queue_[i].size();
        for(int point: queue_[i])
            queue.push_back(point);
    }
}

void beta_BFS(int beta, std::vector<int> &queue, int root){
    /**
     * 顶点索引从0还是1开始：
     * - adja_list从0开始
     * - queue从1开始
    */
    queue.push_back(root);
    if(beta==0) return;

    std::vector<int> pqueue;
    queue.reserve(beta*100);
    pqueue.reserve(beta*100);
    pqueue.push_back(0);

    int layer = 0;
    int begin = 0;
    int last = 1;
    while(layer != beta){
        for(int i = begin; i < last; i++){
            int point = queue[i];
            int parent = pqueue[i];   //parent of point
            for(int search_point: adja_list[point-1]){
                if(parent != search_point+1){
                    queue.push_back(search_point+1);
                    pqueue.push_back(point);
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
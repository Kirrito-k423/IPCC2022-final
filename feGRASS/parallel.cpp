#include "global.h"


/*
使用openMP并行的对vector进行排序
数值相等时选择列号较小的，数值和列号都相等时选择行号较小的
合并时，能被(1<<iteration)整除的进程号，与最近的间隔为1<<(iteration-1)的进程进行合并（1<=iteration<=max_iteration）
进程数可以为任意值，排序后的结果保存在原始vector中
*/
void parallel_sort(vector<vector<double>> &edges) {
    // TODO:尝试使用中间值，使得小于或大于此中间值的数据被分到不同进程，最后直接将结果拼接
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int _num_of_threads = omp_get_num_threads();
        int _total_edges = edges.size();
        int _element_per_thread = (_total_edges + _num_of_threads - 1) / _num_of_threads;
        const int ROW = 0;
        const int COL = 1;
        const int VALUE = 2;
        uint max_iteration = 0;
        int _rank = omp_get_thread_num();
        int start_index = _element_per_thread * _rank;
        int end_index = (_element_per_thread * (_rank + 1) > _total_edges) ? _total_edges: _element_per_thread * (_rank + 1);
        int temp = _num_of_threads;
        while (temp > 1) {
            temp = temp / 2;
            max_iteration++;
        }
        DEBUG_PRINT("rank: %d, edges.size: %ld, start_index: %d, end_index: %d\n", _rank, edges.size(), start_index, end_index);
        //partial_sort(edges.begin(), edges.begin()+end_index, edges.begin()+end_index, compare);
        sort(edges.begin()+start_index, edges.begin()+end_index, compare); // original vector stores sorted result
        DEBUG_PRINT("rank: %d, stable sort finished!\n", _rank);
        for (int iteration = 1; iteration <= max_iteration; iteration++) {
            # pragma omp barrier
            if (_rank % (1 << iteration) == 0) {
                int target_rank = _rank + (1 << (iteration - 1));
                if (target_rank > _num_of_threads) {
                    DEBUG_PRINT("rank: %d break, target rank: %d\n", _rank, target_rank);
                    break;
                }
                vector<double> edge;
                int result_index = start_index;  // 访问原始的vector，修改排序后的值
                end_index = start_index + _element_per_thread * (1 << (iteration-1));
                int target_begin = _element_per_thread * target_rank;
                int target_end = (target_begin+_element_per_thread*(1<<(iteration-1)) >= _total_edges) ? _total_edges : target_begin+_element_per_thread*(1<<(iteration-1));
                int result_end = target_end;
                int index_i = start_index;
                int index_j = target_begin;
                DEBUG_PRINT("rank: %d, target rank: %d, start index: %d, end index: %d, target begin: %d, target end: %d\n", _rank, target_rank, start_index, end_index, target_begin, target_end);
                // 合并两个进程的排序后的结果，保存到原vector中
                vector<vector<double>> copy_of_edges;
                copy_of_edges.assign(edges.begin(), edges.end());
                while(index_i < end_index && index_j < target_end){
                    if(copy_of_edges[index_i][VALUE] > copy_of_edges[index_j][VALUE]  // 数值大
                    // 数值相等时选择列号较小的
                    || ((copy_of_edges[index_i][VALUE] == copy_of_edges[index_j][VALUE]) 
                        && (copy_of_edges[index_i][COL] < copy_of_edges[index_j][COL]))
                    // 数值和列号都相等时选择行号较小的
                    || ((copy_of_edges[index_i][VALUE] == copy_of_edges[index_j][VALUE]) 
                        && (copy_of_edges[index_i][COL] == copy_of_edges[index_j][COL]) 
                        && (copy_of_edges[index_i][ROW] < copy_of_edges[index_j][ROW]))
                    ){
                        // 因为是不同的点，行号和列号不会同时相等
                        // 这里使用memory copy并不明显的比直接赋值快
                        memcpy(edges[result_index].data(), copy_of_edges[index_i].data(), sizeof(double) * 4);
                        index_i++;
                        result_index++;
                    }
                    else{
                        memcpy(edges[result_index].data(), copy_of_edges[index_j].data(), sizeof(double) * 4);
                        index_j++;
                        result_index++;
                    }
                }
                while(index_i < end_index){
                    memcpy(edges[result_index].data(), copy_of_edges[index_i].data(), sizeof(double) * 4);
                    index_i++;
                    result_index++;
                }
                while(index_j < target_end){
                    memcpy(edges[result_index].data(), copy_of_edges[index_j].data(), sizeof(double) * 4);
                    index_j++;
                    result_index++;
                }
                copy_of_edges.clear();
            }
        }
    }
}
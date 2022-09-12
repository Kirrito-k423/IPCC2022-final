#include "global.h"


/*
使用openMP并行的对vector进行排序
数值相等时选择列号较小的，数值和列号都相等时选择行号较小的
合并时，能被(1<<iteration)整除的进程号，与最近的间隔为1<<(iteration-1)的进程进行合并（1<=iteration<=max_iteration）
进程数可以为任意值，排序后的结果保存在原始vector中
*/
void parallel_sort(vector<vector<double>> &edges) {
    // TODO:尝试使用中间值，使得小于或大于此中间值的数据被分到不同进程，最后直接将结果拼接
    int _num_of_threads = omp_get_num_threads();
    //printf("using thread num: %d", _num_of_threads);
    int _total_edges = edges.size();
    int _element_per_thread = (_total_edges + _num_of_threads - 1) / _num_of_threads;
    const int ROW = 0;
    const int COL = 1;
    const int VALUE = 2;
    #pragma omp parallel 
    {
        uint max_iteration = 0;
        int _rank = omp_get_thread_num();
        //printf("rank: %d\n", _rank);
        int start_index = _element_per_thread * _rank;
        int end_index = (_rank == _num_of_threads) ? _total_edges : _element_per_thread * (_rank + 1);
        int temp = _num_of_threads;
        while (temp > 1) {
            temp = temp / 2;
            max_iteration++;
        }
        #pragma omp for
        for (int i = 0; i < _num_of_threads; i++) {
            sort(edges.begin()+start_index, edges.begin()+end_index, compare2); // original vector stores sorted result
            for (int iteration = 0; iteration < max_iteration; iteration++) {
                if (_rank % (1 << iteration) == 0) {
                    int target_rank = _rank + (1 << iteration + 1 << (iteration - 1));
                    if (target_rank > _num_of_threads) {
                        break;
                    }
                    vector<double> edge;
                    int result_index = start_index;  // 访问原始的vector，修改排序后的值
                    end_index = start_index + _element_per_thread * (1 << iteration);
                    int target_begin = _element_per_thread * target_rank;
                    int target_end = (target_rank + (1 << (iteration-1)) >= _num_of_threads) ? _total_edges : target_begin + _element_per_thread * (1 << iteration);
                    int result_end = target_end;
                    int index_i = start_index;
                    int index_j = target_begin;
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
                        memcpy(edges[result_index].data(), copy_of_edges[index_i].data(), sizeof(double) * 4);
                        index_j++;
                        result_index++;
                    }
                    copy_of_edges.clear();
                }
            }
        }
    }
}
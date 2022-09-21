#include "global.h"
int cmp_by_effw(const void *a, const void *b) {
    return (*(edge_t *)a) < (*(edge_t *)b);
}

void psrs(vector<edge_t> &arr, int p, comparison_fn_t cmp) {
    //struct timeval startTime, endTime;
    //gettimeofday(&startTime, NULL);
    int size = arr.size();
    edge_t **local_arr = (edge_t **)malloc(p * sizeof(edge_t *));      // 局部数组
    edge_t **sample_arr = (edge_t **)malloc(p * sizeof(edge_t *));     // 每个进程的采样结果
    edge_t *pivots = (edge_t *)malloc((p - 1) * sizeof(edge_t));       // 分区主元
    int **index_arr = (int **)malloc(p * sizeof(int *));               // 按照主元划分之后，每个进程中，每块的下标索引
    int **size_arr = (int **)malloc(p * sizeof(int *));                // 按照主元划分之后，每个进程中，每块的大小
    int *sorted_size_arr = (int *)malloc(p * sizeof(int));             // 全局交换后，每个进程保存的元素数目
    edge_t **sorted_ptr_arr = (edge_t **)malloc(p * sizeof(edge_t *)); // 全局交换之后，每个进程的数组索引
    //printTime("psrs sort:before openMP");
    #pragma omp parallel num_threads(p)
    {
        // divide to each processor
        int tid = omp_get_thread_num();
        int ele_per_processor = (size + p - 1) / p;
        int start = tid * ele_per_processor;
        int end = (tid + 1) * ele_per_processor > size ? size : (tid + 1) * ele_per_processor;
        int length = end - start;
        local_arr[tid] = (edge_t *)malloc(length * sizeof(edge_t));
        // TODO: local_arr可以被优化掉
        memcpy(local_arr[tid], arr.data() + start, length * sizeof(edge_t));

        qsort(local_arr[tid], length, sizeof(edge_t), cmp);
        // store samples
        sample_arr[tid] = (edge_t *)malloc(p * sizeof(edge_t));
        for (int i = 0; i < p; i++) {
            sample_arr[tid][i] = local_arr[tid][i * length / p];
            DEBUG_PRINT("tid: %d, sorted index: %d, eff_w: %lf\n", tid, i * length / p, local_arr[tid][i * length / p].eff_w);
        }
        #pragma omp barrier
        if (tid == 0) {
            // main processor merge sort samples
            edge_t *sorted_samples = (edge_t *)malloc(p * p * sizeof(edge_t));
            for (int i = 0; i < p; i++) {
                memcpy((void *)(sorted_samples + p * i), (void *)sample_arr[i], p * sizeof(edge_t));
            }
            qsort(sorted_samples, p * p, sizeof(edge_t), cmp);
            for (int i = 0; i < p - 1; i++) {
                pivots[i] = sorted_samples[(i + 1) * p];
                DEBUG_PRINT("sorted_samples[%d]: %d eff_w: %lf\n", (i + 1) * p, pivots[i], pivots[i].eff_w);
            }
        }
        index_arr[tid] = (int *)malloc(sizeof(int) * p); // 划分的下标索引
        memset(index_arr[tid], 0, sizeof(int) * p);
        size_arr[tid] = (int *)malloc(sizeof(int) * p); // 划分的大小
        memset(size_arr[tid], 0, sizeof(int) * p);
        int index_i = 0; // 遍历局部数组
        int index_j = 0; // 遍历主元
        int count = 0;
        #pragma omp barrier
        // 每个进程根据pivots将自己进程内分配的数据进行划分
        while (index_i < length && index_j < p - 1) {
            // 目的是降序排列，传入的cmp函数为参数1<参数2时为true（对应快排，满足时交换位置）
            // 这里应该是pivot[index_j]对应的值小于local_arr_[index_i]且大于local_arr_[index_i+1]
            if (cmp(&local_arr[tid][index_i], &pivots[index_j])) {
                // pivot[index_j] > local_arr[tid][index_i]
                DEBUG_PRINT("tid %d: pivots[%d]: %lf > local_arr[%d][%d]: %lf, count: %d\n", tid, index_j, pivots[index_j].eff_w, tid, index_i, local_arr[tid][index_i].eff_w, count);
                size_arr[tid][index_j] = count;
                count = 0;
                index_j++;
                index_arr[tid][index_j] = index_i;
            } else if (cmp(&pivots[index_j], &local_arr[tid][index_i])) {
                // local_arr[index_i] == pivot[index_j]  or local_arr[index_i] > pivot[index_j]
                index_i++;
                count++;
            } else if (local_arr[tid][index_i] == pivots[index_j]) {
                DEBUG_PRINT("local_arr[%d][%d]] %lf == pivots[%d] %lf\n", tid, index_i, local_arr[tid][index_i].eff_w, index_j, pivots[index_j].eff_w);
                index_i++;
                count++;
            }
        }
        if (index_i != length - 1) {
            // 对应pivots已经遍历完，但是local_arr[tid]没有遍历完，index_j一定为p-1
            // 同时pivots中最后一个元素，大于local_arr[tid]中剩余所有元素
            // local_arr[tid]中剩余的一部分，被主元划分到最后一个分区
            DEBUG_PRINT("tid: %d, remain %d of %d, index_i: %d, index_j: %d start: %d, end: %d\n", tid, length - index_i, length, index_i, index_j, start, end);
            size_arr[tid][index_j] = length - index_i;
        }
        if (index_j != p - 1) {
            // TODO: case 4 未经历此分支
            // 说明pivots没有遍历完，即pivots中存在一个元素，小于local_arr[tid]中所有元素
            DEBUG_PRINT("pivots left, tid: %d, index_i: %d, index_j: %d length: %d, start: %d, end: %d\n", tid, index_i, length, index_j, start, end);
            for (int i = index_j; i < p - 1; i++) {
                index_arr[tid][i] = index_i;
                size_arr[tid][i] = 0;
            }
        }
        // int tmp = 0;
        // for (int i = 0; i < p; i++) {
        //     tmp += size_arr[tid][i];
        //     DEBUG_PRINT("tid %d[%d]: index: %d, size: %d, length: %d\n", tid, i, index_arr[tid][i], size_arr[tid][i], length);
        // }
        free(sample_arr[tid]);
        // TODO:使用优先队列进行多路合并（堆排序）
        int sorted_size = 0;
        #pragma omp barrier
        for (int i = 0; i < p; i++) {
            // 统计应该分配到自己进程的元素数目
            sorted_size += size_arr[i][tid];
        }
        sorted_size_arr[tid] = sorted_size;
        DEBUG_PRINT("tid: %d sorted size: %d\n", tid, sorted_size);
        sorted_ptr_arr[tid] = (edge_t *)malloc(sizeof(edge_t) * sorted_size);
        int ele_count = 0; // 将各个进程按照主元划分之后的分区中属于自己进程的部分，拷贝到自己的变量中
        for (int i = 0; i < p; i++) {
            DEBUG_PRINT("tid: %d, sorted_ptr_arr[%d] + %d, local_arr[%d] + index_arr[%d][%d](%d), size_arr[%d][%d](%d)\n", tid, tid, ele_count, i, i, tid, index_arr[i][tid], i, tid, size_arr[i][tid]);
            memcpy(sorted_ptr_arr[tid] + ele_count, local_arr[i] + index_arr[i][tid], size_arr[i][tid] * sizeof(edge_t));
            ele_count += size_arr[i][tid];
        }
        qsort(sorted_ptr_arr[tid], sorted_size, sizeof(edge_t), cmp);
        int pos = 0; // 此进程排序后的数组，在全局数组进行拼接时的开始位置
        for (int i = 0; i < tid; i++) {
            pos += sorted_size_arr[i];
        }
        memcpy(arr.data() + pos, sorted_ptr_arr[tid], sizeof(edge_t) * sorted_size);
    }
    for (int i = 0; i < p; i++) {
    // 这里不能在openMP内free（除非加barrier），因为memory copy时会有冲突
        free(size_arr[i]);
        free(index_arr[i]);
        free(local_arr[i]);
        free(sorted_ptr_arr[i]);
    }
    free(local_arr);
    free(sample_arr);
    free(index_arr);
    free(size_arr);
    free(sorted_size_arr);
    free(sorted_ptr_arr);
}

// template <typename T>
// void psrs(vector<T> &arr, int p, comparison_fn_t cmp) {
//     int size = arr.size();
//     T **local_arr = (T **)malloc(p * sizeof(T *));         // 局部数组
//     T **sample_arr = (T **)malloc(p * sizeof(T *));        // 每个进程的采样结果
//     T *pivots = (T *)malloc((p - 1) * sizeof(T));          // 分区主元
//     int **index_arr = (int **)malloc(p * sizeof(int *));   // 按照主元划分之后，每个进程中，每块的下标索引
//     int **size_arr = (int **)malloc(p * sizeof(int *));    // 按照主元划分之后，每个进程中，每块的大小
//     int *sorted_size_arr = (int *)malloc(p * sizeof(int)); // 全局交换后，每个进程保存的元素数目
//     T **sorted_ptr_arr = (T **)malloc(p * sizeof(T *));    // 全局交换之后，每个进程的数组索引
// #pragma omp prallel num_threads(p)
//     {
//         // divide to each processor
//         int tid = omp_get_thread_num();
//         int ele_per_processor = (size + p - 1) / p;
//         int start = tid * ele_per_processor;
//         int end = (tid + 1) * ele_per_processor > size ? size : (tid + 1) * ele_per_processor;
//         int length = start - end;
//         local_arr[tid] = (T *)malloc(length * sizeof(T));
//         // TODO: local_arr可以被优化掉
//         memcpy(local_arr[tid], arr.data() + start, length * sizeof(T));
//         qsort(local_arr[tid], length, sizeof(T), cmp);
//         // store samples
//         sample_arr[tid] = (T *)malloc(p * sizeof(T));
//         for (int i = 0; i < p; i++) {
//             sample_arr[tid][i] = local_arr[tid][i * length / p];
//         }
// #pragma omp barrier
//         if (tid == 0) {
//             // main processor merge sort samples
//             T *sorted_samples = (T *)malloc(p * p * sizeof(T));
//             for (int i = 0; i < p; i++) {
//                 memcpy((void*)sorted_samples[p * i], (void*)sample_arr[i], p * sizeof(T));
//             }
//             qsort(sorted_samples, p * p, sizeof(T), cmp);
//             for (int i = 0; i < p - 1; i++) {
//                 pivots[i] = sorted_samples[(i + 1) * p];
//             }
//         }
// #pragma omp barrier
//         index_arr[tid] = (int *)malloc(sizeof(int) * p); // 划分的下标索引
//         memset(index_arr[tid], 0, sizeof(int) * p);
//         size_arr[tid] = (int *)malloc(sizeof(int) * p); // 划分的大小
//         memset(size_arr[tid], 0, sizeof(int) * p);
//         int index_i = 0; // 遍历局部数组
//         int index_j = 0; // 遍历主元
//         int count = 0;
//         int partition = 0; // 按主元划分之后的分块索引，从0开始
//         // 每个进程根据pivots将自己进程内分配的数据进行划分
//         while (index_i < length && index_j < p ) {
//             // 目的是降序排列，传入的cmp函数为参数1<参数2时为true（对应快排，满足时交换位置）
//             // 这里应该是pivot[index_j]对应的值小于local_arr_[index_i]且大于local_arr_[index_i+1]
//             // cmp函数应该能够保证不会出现相等的情况
//             if (cmp(&local_arr[index_i], &pivots[index_j])) {
//                 // pivot[index_j] > local_arr[index_i]
//                 size_arr[tid][partition] = count;
//                 count = 0;
//                 partition++; index_j++;
//                 if (partition == p) {
//                     break;
//                 }
//                 index_arr[tid][partition] = index_i;
//             } else if (local_arr[index_i] == pivots[index_j] || cmp(&pivots[index_j], &local_arr[index_i])) {
//                 // pivot[index_j] == local_arr[index_i] or pivot[index_j] < local_arr[index_i]
//                 index_i++;
//                 count++;
//             }
//         }
//         if(index_i != length - 1){
//             // local_arr中剩余的一部分，被主元划分到最后一个分区
//             size_arr[tid][partition] = length - index_i;
//         }
//         DEBUG_PRINT("tid: %d, ")
//         free(local_arr[tid]);
//         free(sample_arr[tid]);
//         // TODO:使用优先队列进行多路合并（堆排序）
//         int sorted_size = 0;
//         for (int i = 0; i < p; i++) {
//             sorted_size += size_arr[i][tid];
//         }
//         sorted_size_arr[tid] = sorted_size;
//         sorted_ptr_arr[tid] = (T *)malloc(sizeof(T) * sorted_size);
//         int ele_count = 0; // 将各个进程按照主元划分之后的分区中属于自己进程的部分，拷贝到自己的变量中
//         for (int i = 0; i < p; i++) {
//             memcpy((void*)(sorted_ptr_arr[tid] + ele_count), (const void *)local_arr[tid][index_arr[i][tid], size_arr[i][tid]*sizeof(T));
//             ele_count += sorted_size_arr[i];
//         }
//         free(index_arr[tid]);
//         free(size_arr[tid]);
//         qsort(sorted_ptr_arr[tid], ele_count, sizeof(T), cmp);
//         int pos = 0; // 此进程排序后的数组，在全局数组进行拼接时的开始位置
//         for (int i = 0; i < tid; i++) {
//             pos += sorted_size_arr[i];
//         }
//         memcpy(arr.data() + pos, sorted_ptr_arr, sizeof(T) * sorted_size);
//         free(sorted_ptr_arr[tid]);
//     }
//     free(local_arr);
//     free(sample_arr);
//     free(index_arr);
//     free(size_arr);
//     free(sorted_size_arr);
//     free(sorted_ptr_arr);
// }
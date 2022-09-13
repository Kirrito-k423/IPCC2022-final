template <typename T>
T *merge(T *arr1, int len1, T *arr2, int len2, comparison_fn_t cmp) {
    T *arr = (T *)malloc((len1 + len2) * sizeof(T));
    int i, j, idx;
    idx = i = j = 0;
    for (; i < len1 && j < len2;) {
        // if (arr[p2 + j] > arr[p1 + i]) {     //降序排序，相同时，优先arr[p1 + i]（索引小）
        if (cmp(&arr1[i], &arr2[j])) {
            arr[idx++] = arr2[j];
            j++;
        } else {
            arr[idx++] = arr1[i];
            i++;
        }
    }
    if (i >= len1) {
        for (; j < len2; j++) {
            arr[idx++] = arr2[j];
        }
    } else {
        for (; i < len1; i++) {
            arr[idx++] = arr1[i];
        }
    }
    free(arr1);
    free(arr2);
    return arr;
}

/**
 * p: 线程数目，必须为2的幂
 */
template <typename T>
void p_mergesort(vector<T> &arr, int p, comparison_fn_t cmp) {
    // T *p_mergesort(vector<T> &arr, int p, bool (*cmp)(const T &a, const T &b)) {
    int n = arr.size();
    int blk_size[p];
    int offset[p + 1];
    memset(offset, 0, sizeof(offset));
    offset[p] = n;
    T *arr_[p];
    for (int i = 0; i < p; i++) {
        blk_size[i] = n / p;
        // if(tid==p-1){    //(p-1)*q + q+r
        //     blk_size[tid] += n % p;
        // }
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
        // printf("%d: %d %d\n", tid, offset[tid], blk_size[tid]);

        arr_[tid] = (T *)malloc(blk_size[tid] * sizeof(T));
        memcpy(arr_[tid], arr.data() + offset[tid], blk_size[tid] * sizeof(T));

        // sort
        qsort(arr_[tid], blk_size[tid], sizeof(T), cmp);
        // printf("%d: sort finish\n", tid);

#pragma omp barrier
        // binary tree merge
        int i = 1;
        while (i < p) {         // log p轮
            int pair = tid ^ i; // tid需要和pair合并
            i <<= 1;
            if (tid % i == 0) {
                // printf("(%d, %d) round %d\n", tid, pair, i>>1);
                arr_[tid] = merge(arr_[tid], blk_size[tid], arr_[pair], blk_size[pair], cmp); //合并的都是相邻的
                blk_size[tid] += blk_size[pair];
            }
#pragma omp barrier
        }
    }
    memcpy(arr.data(), arr_[0], n * sizeof(T));
}
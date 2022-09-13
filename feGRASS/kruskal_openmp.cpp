# include "global.h"


edge * merge_pair_msf(int vertices, edge *local_edges, int len1, edge *pair_edges, int len2, int *edge_cnt, comparison_fn_t cmp){
    DSU dsu = DSU(vertices);
    edge *msf_result = (edge *) malloc((len1 + len2) * sizeof(edge));
    /* 
    Local_edges and pair_edges are already sorted in order, 
    so we can check each edge in these two array one by one,
    to check if they form a loop
    */
    int i, j;
    i = j = 0;
    while(i < len1 && j < len2){
        if(cmp(&local_edges[i], &pair_edges[j])){
            int u = int(local_edges[i].u);
            int v = int(local_edges[i].v);
            if(!dsu.same_set(u, v)){
                msf_result[*edge_cnt++] = local_edges[i];
                dsu.unite(u, v);
                *edge_cnt++;
            }
            i++;
        } else{
            int u = int(pair_edges[j].u);
            int v = int(pair_edges[j].v);
            if(!dsu.same_set(u, v)){
                msf_result[*edge_cnt++] = local_edges[i];
                dsu.unite(u, v);
                *edge_cnt++;
            }
            j++;
        }
    }
    while(i < len1){
        int u = int(local_edges[i].u);
        int v = int(local_edges[i].v);
        if(!dsu.same_set(u, v)){
            msf_result[*edge_cnt++] = local_edges[i];
            dsu.unite(u, v);
            *edge_cnt++;
        }
        i++;
    }
    while(j < len2){
        int u = int(pair_edges[j].u);
        int v = int(pair_edges[j].v);
        if(!dsu.same_set(u, v)){
            msf_result[*edge_cnt++] = local_edges[i];
            dsu.unite(u, v);
            *edge_cnt++;
        }
        j++;
    }
    return msf_result;
}



void p_kruskal(int vertices, vector<edge> &edge_matrix, vector<edge_t> &spanning_tree, int p){
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);
    // parallel portion using openMP

    /*
    1. load edges distributed to each processor
    2. every processor find local minimum spanning forest
    3. merge from each pair of minimum spanning forest
    */
    int n = edge_matrix.size();
    // store edge number of spanning forest for each processor
    int *msf_edges_num = (int *) malloc(p * sizeof(int));
    memset(msf_edges_num, 0, p * sizeof(int));
    // store address of spanning forest for each processor
    edge** msf_edges_ptr = (edge**)malloc(p * sizeof(edge*));
    #pragma opm parallel num_threads(p)
    {
        int tid = omp_get_thread_num();
        int vertices_per_processor = (vertices + p - 1) / p;
        int start_vertices_index = tid * vertices_per_processor;
        int end_vertices_index = (tid + 1) * vertices_per_processor > n ? n : (tid + 1) * vertices_per_processor;
        edge *local_edges = (edge *) malloc(n * sizeof(edge));
        int used_edges = 0;
        for(int i = 0; i < n; i++){
            if((edge_matrix[i].u >= start_vertices_index && edge_matrix[i].u < end_vertices_index)
            || (edge_matrix[i].v >= start_vertices_index && edge_matrix[i].v < end_vertices_index))
            local_edges[used_edges++] = edge_matrix[i];
        }
        DSU dsu = DSU(n);
        qsort(local_edges, used_edges, sizeof(edge), cmp);
        edge *local_spanning_tree = (edge *) malloc(used_edges * sizeof(edge));
        int sf_edge_cnt = 0;
        for (int i=0; i<used_edges; i++) {
            int u = int(local_edges[i].u);
            int v = int(local_edges[i].v);
            if(!dsu.same_set(u, v)){
                local_spanning_tree[sf_edge_cnt++] = local_edges[i];
                dsu.unite(u, v);
                sf_edge_cnt++;
            }
        }
        free(local_edges);
        msf_edges_num[tid] = sf_edge_cnt;
        msf_edges_ptr[tid] = local_spanning_tree;
        int i = 0;
        while(i < p) {
            #pragma omp barrier
            int pair = tid ^ i;
            i <<= 1;
            if (tid % i == 0){
                // merge
                int *merged_edge_cnt = 0;
                merge_pair_msf(vertices, msf_edges_ptr[tid], msf_edges_num[tid], msf_edges_ptr[pair], msf_edges_num[pair], msf_edges_num, cmp);
                msf_edges_num[tid] = *merged_edge_cnt;
                free(msf_edges_ptr[pair]);
            }
        }
    }
    memcpy(spanning_tree.data(), msf_edges_ptr, (vertices-1)*sizeof(edge));
    printTime("construct spanning tree using openMP");
}
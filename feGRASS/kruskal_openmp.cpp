#include "global.h"

int cmp(const void *a, const void *b) {
    return ((edge_t *)a)->eff_w < ((edge_t *)b)->eff_w;
}

//disjoint set union
class DSU
{
private:
    int size;
    int *parent;
    int *rank;
public:
    DSU(int size): size(size){
        parent = new int[size];
        rank = new int[size];
        for(int i=0; i<size; i++){
            parent[i] = i;
            rank[i] = 0;
        }
    }
    ~DSU(){
        delete []parent;
        delete []rank;
    }
    int find_root(int id) {
        while (id != parent[id]) {
            id = parent[id];
        }

        return id;
    }

    bool same_set(int id1, int id2) {
        return find_root(id1) == find_root(id2);
    }

    void unite(int id1, int id2) {
        id1 = find_root(id1);
        id2 = find_root(id2);

        if (rank[id1] < rank[id2]) std::swap(id1, id2);

        parent[id2] = id1;
        if (rank[id1] == rank[id2]) ++rank[id1];
    }
};


edge * merge_pair_msf(int vertices, edge *local_edges, int len1, edge *pair_edges, int len2, int *edge_cnt, comparison_fn_t cmp){
    DSU dsu = DSU(vertices+1);
    edge *msf_result = (edge *) malloc((len1 + len2) * sizeof(edge));
    /* 
    Local_edges and pair_edges are already sorted in order, 
    so we can check each edge in these two array one by one,
    to check if they form a loop
    */
    int i, j;
    i = j = 0;
    while(i < len1 && j < len2){
        if(local_edges[i].eff_w == pair_edges[j].eff_w){
            if((local_edges[i].v < pair_edges[j].v)
            || (local_edges[i].v == pair_edges[j].v) && local_edges[i].u < pair_edges[j].u){
                int u = int(local_edges[i].u);
                int v = int(local_edges[i].v);
                if(!dsu.same_set(u, v)){
                    msf_result[*edge_cnt] = local_edges[i];
                    dsu.unite(u, v);
                    (*edge_cnt) += 1;
                }
                i++;
            } else{
                int u = int(pair_edges[j].u);
                int v = int(pair_edges[j].v);
                if(!dsu.same_set(u, v)){
                    msf_result[*edge_cnt] = pair_edges[j];
                    dsu.unite(u, v);
                    (*edge_cnt) += 1;
                }
                j++;
            }
            continue;
        }
        if(!cmp(&local_edges[i], &pair_edges[j])){
            int u = int(local_edges[i].u);
            int v = int(local_edges[i].v);
            if(!dsu.same_set(u, v)){
                msf_result[*edge_cnt] = local_edges[i];
                dsu.unite(u, v);
                (*edge_cnt) += 1;
            }
            i++;
        } else{
            int u = int(pair_edges[j].u);
            int v = int(pair_edges[j].v);
            if(!dsu.same_set(u, v)){
                msf_result[*edge_cnt] = pair_edges[j];
                dsu.unite(u, v);
                (*edge_cnt) += 1;
            }
            j++;
        }
    }
    while(i < len1){
        int u = int(local_edges[i].u);
        int v = int(local_edges[i].v);
        if(!dsu.same_set(u, v)){
            msf_result[*edge_cnt] = local_edges[i];
            dsu.unite(u, v);
            (*edge_cnt) += 1;
        }
        i++;
    }
    while(j < len2){
        int u = int(pair_edges[j].u);
        int v = int(pair_edges[j].v);
        if(!dsu.same_set(u, v)){
            msf_result[*edge_cnt] = pair_edges[j];
            dsu.unite(u, v);
            (*edge_cnt) += 1;
        }
        j++;
    }
    return msf_result;
}



void p_kruskal(int vertices, vector<edge_t> &edge_matrix, vector<edge_t> &spanning_tree, int p){
    //struct timeval startTime, endTime;
    //gettimeofday(&startTime, NULL);
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
    edge_t** msf_edges_ptr = (edge_t**)malloc(p * sizeof(edge_t*));
    #pragma omp parallel num_threads(p)
    {
        int tid = omp_get_thread_num();
        int vertices_per_processor = (vertices + p - 1) / p;
        int start_vertices_index = tid * vertices_per_processor;
        int end_vertices_index = (tid + 1) * vertices_per_processor > n ? n : (tid + 1) * vertices_per_processor;
        edge_t *local_edges = (edge_t *) malloc(n * sizeof(edge_t));
        int used_edges = 0;
        for(int i = 0; i < n; i++){
            if((edge_matrix[i].u >= start_vertices_index && edge_matrix[i].u < end_vertices_index)
            || (edge_matrix[i].v >= start_vertices_index && edge_matrix[i].v < end_vertices_index))
            local_edges[used_edges++] = edge_matrix[i];
        }
        //DEBUG_PRINT("tid: %d get local edges\n", tid);
        DSU dsu = DSU(n+1);
        qsort(local_edges, used_edges, sizeof(edge_t), cmp);
        edge_t *local_spanning_tree = (edge_t *) malloc(used_edges * sizeof(edge_t));
        int sf_edge_cnt = 0;
        for (int i=0; i<used_edges; i++) {
            int u = int(local_edges[i].u);
            int v = int(local_edges[i].v);
            if(!dsu.same_set(u, v)){
                local_spanning_tree[sf_edge_cnt++] = local_edges[i];
                dsu.unite(u, v);
            }
        }
        //DEBUG_PRINT("tid: %d find local spanning forest\n", tid);
        free(local_edges);
        msf_edges_num[tid] = sf_edge_cnt;
        msf_edges_ptr[tid] = local_spanning_tree;
        int i = 1;
        while(i < p) {
            #pragma omp barrier
            int pair = tid ^ i;
            i <<= 1;
            if (tid % i == 0){
                // merge
                //DEBUG_PRINT("tid: %d, pair tid: %d, start to merge\n", tid, pair);
                int merged_edge_cnt = 0;
                msf_edges_ptr[tid] = merge_pair_msf(vertices, msf_edges_ptr[tid], msf_edges_num[tid], msf_edges_ptr[pair], msf_edges_num[pair], &merged_edge_cnt, cmp);
                msf_edges_num[tid] = merged_edge_cnt;
                free(msf_edges_ptr[pair]);
            }
        }
    }
    DEBUG_PRINT("vertices: %d, edge num: %d\n", vertices, msf_edges_num[0]);
    memcpy(spanning_tree.data(), msf_edges_ptr[0], (vertices-1)*sizeof(edge_t));
    //printTime("construct spanning tree using openMP");
}
#include "global.h"

vector<vector<int>> adja_list; //给之后的BFS使用。adja_list_w只在DFS_traversal中使用
vector<vector<node_t>> adja_list_w;
double *dis;
int *parent;
int *no_weight_dis;

int cmp_by_index(const void *a, const void *b) {
    int u_less = ((edge_t *)a)->u < ((edge_t *)b)->u;
    int u_equal = ((edge_t *)a)->u == ((edge_t *)b)->u;
    int v_less = ((edge_t *)a)->v < ((edge_t *)b)->v;
    return u_equal ? v_less : u_less;
}


/**
 * 根据树的边集表示，创建邻接表
*/
void create_adja_list(vector<edge_t> &tree, vector<vector<node_t>> &adja_list_w){
    int edge_num = tree.size();
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);
    int tree_size = tree.size();
    tree.resize(2 * tree_size);
    // duplicate edges, the back half edges aim to reverted direction
    memcpy(tree.data() + tree_size, tree.data(), tree_size * sizeof(edge_t));
    int p = CREATE_ADJA_THREADS;
    #pragma omp parallel num_threads(p)
    {
        int tid = omp_get_thread_num();
        // reverted direction of the back half edges parallelly
        int edge_start = edge_num + tid * edge_num / p;
        int edge_end = edge_num + (tid + 1) * edge_num / p;
        // DEBUG_PRINT("tid: %d, revert edges start: %d, end: %d\n", tid, edge_start, edge_end);
        for(int index = edge_start; index < edge_end; index++){
            // revert direction
            int tmp = tree[index].u;
            tree[index].u = tree[index].v;
            tree[index].v = tmp;
        }
        // traverse the whole edges (with both direction), and append edges to assigned nodes
        int node_start = tid * M / p;
        int node_end = (tid+1) * M / p;
        // DEBUG_PRINT("tid: %d, construct ajda_list, node start: %d, end: %d\n", tid, node_start, node_end);
        #pragma omp barrier
        node_t node;    //[point_idx, edge_weight]
        for(int i=0; i<2*edge_num; i++){
            int node_index = tree[i].u-1;
            if(node_index < node_end && node_index >= node_start){
                node.u = tree[i].v - 1;
                node.w = tree[i].w;
                adja_list_w[node_index].push_back(node);
                adja_list[node_index].push_back(tree[i].v - 1);
            }
        }
    }
    tree.resize(tree_size);
}

/**
 * DFS遍历生成树获得：
 * 1. 每个点到根节点的距离
 * 2. 每个点的parent数组
 * 3. 每个点到根节点的无权重距离
 * 根节点使用largest_volume_point
*/
void DFS_traversal(vector<vector<node_t>> &adja_list_w, double *dis, int *parent, int *no_weight_dis){
    //init distance array
    memset(dis, 0, sizeof(double)*M);
    //init parent array
    memset(parent, 0, sizeof(int)*M);
    parent[largest_volume_point-1] = largest_volume_point-1;
    //no weight distance
    memset(no_weight_dis, -1, sizeof(int)*M);
    no_weight_dis[largest_volume_point-1] = 0;
    //DFS
    stack<int> process;
    int point;
    process.push(largest_volume_point-1);
    while(!process.empty()){
        point=process.top();
        process.pop();
        for(int i=0; i<adja_list_w[point].size(); i++){
            int search_point = adja_list_w[point][i].u;
            if(no_weight_dis[search_point]==-1){   //first search
                process.push(search_point);
                dis[search_point] = dis[point] + 1.0/adja_list_w[point][i].w;    //add edge weight
                parent[search_point] = point;
                no_weight_dis[search_point] = no_weight_dis[point] + 1;
            }
        }
    }
}

void debug_print(double *dis, int *parent, int *no_weight_dis){
    printf("largest_volume_point: %d\n", largest_volume_point-1);
    printf("%4s %4s %4s %4s\n", "i", "pare", "dis", "wdis");
    for(int i=0; i<M; i++){
        printf("%4d %4d %4d %4.1f\n", i, parent[i], no_weight_dis[i], dis[i]);
    }
}

void debug_print_adja(vector<vector<node_t>> &adja_list){
    printf("adja: \n");
    for(int i=0; i<M; i++){
        printf("%2d: ", i);
        for(int j=0; j<adja_list[i].size(); j++){
            printf("%d ", int(adja_list[i][j].u));
        }
        printf("\n");
    }
}

/**
 * 打印i, j到LCA上途径的点
*/
void debug_print_path(int i, int j, int *parent, int *no_weight_dis){
    if(no_weight_dis[i] < no_weight_dis[j]){
        int tmp = i;
        i = j;
        j = tmp;
    }
    //keep dis[i] >= dis[j]
    printf("%4d %4d\n", i, j);
    int delta = no_weight_dis[i] - no_weight_dis[j];
    while(delta--){
        i = parent[i];
        printf("%4d\n", i);
    }

    if(i == j) return;

    while(parent[i]!=parent[j]){
        i = parent[i];
        j = parent[j];
        printf("%4d %4d\n", i, j);
    }
    printf("  %4d\n", parent[i]);
}

/**
 * 计算i, j的最近公共祖先节点
 * 使用parent数组
*/
int get_LCA(int i, int j, int *parent, int *no_weight_dis){
    if(no_weight_dis[i] < no_weight_dis[j]){
        int tmp = i;
        i = j;
        j = tmp;
    }
    //keep dis[i] >= dis[j]
    int delta = no_weight_dis[i] - no_weight_dis[j];
    while(delta--){
        i = parent[i];
    }

    if(i == j) return i;    //如果i和j位于一条路径上

    while(parent[i]!=parent[j]){
        i = parent[i];
        j = parent[j];
    }
    return parent[i];
}

/**
 * 计算每个非树边在生成树上的等效电阻
 * 结果存在copy_off_tree_edge
 * 
 * 格式：
 * spanning_tree: [point_index1, point_index2, W_eff, W_ij]
 * off_tree_edge: subset of spanning_tree
 * copy_off_tree_edge: [point_index1, point_index2, R_ij, W_ij]
*/
void caculate_resistance(vector<edge_t> &spanning_tree, vector<edge_t> &off_tree_edge, vector<edge_t> &copy_off_tree_edge){
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);

    adja_list.resize(M);
    adja_list_w.resize(M);
    dis = (double *)malloc(M*sizeof(double));   //每个点到根节点距离
    parent = (int *)malloc(M*sizeof(int));         //每个点父节点数组
    no_weight_dis = (int *)malloc(M*sizeof(int));  //每个点到根节点无权重距离

    //创建生成树的邻接表
    create_adja_list(spanning_tree, adja_list_w);
    printTime("resistance: create adja list")
    // debug_print_adja(adja_list);
    
    //DFS遍历，得到dis,parent等
    DFS_traversal(adja_list_w, dis, parent, no_weight_dis);
    printTime("resistance: DFS traverse spanning tree")
    // debug_print(dis, parent, no_weight_dis);
    
    // printf("path: \n");
    // debug_print_path(1859, 3044, parent, no_weight_dis);
    
    int n = off_tree_edge.size();
    const int p = NUM_THREADS;
    int blk_size[p];
    int offset[p + 1];
    memset(offset, 0, sizeof(offset));
    offset[p] = n;
    for (int i = 0; i < p; i++) {
        blk_size[i] = n / p;
        if (i < n % p) { // r*(q+1) + (p-r)*q
            blk_size[i] += 1;
        }
    }
    copy_off_tree_edge.resize(n);
    #pragma omp parallel num_threads(p)
    {
        const int tid=omp_get_thread_num();

        for (int i = 0; i < tid; i++) {
            offset[tid] += blk_size[i];
        }
        for (int i=offset[tid]; i<offset[tid] + blk_size[tid]; i++) {
            edge_t edge = off_tree_edge[i];
            int edge_point1 = edge.u-1;
            int edge_point2 = edge.v-1;
            //树上计算等效电阻简单方式
            int LCA_point = get_LCA(edge_point1, edge_point2, parent, no_weight_dis);
            // printf("LCA(%d, %d)=%d\n", edge_point1, edge_point2, LCA_point);
            // int d_ = no_weight_dis[edge_point1] + no_weight_dis[edge_point2] - 2*no_weight_dis[LCA_point];
            // printf("no_weight_dis(%d, %d)=%d\n", edge_point1, edge_point2, d_);
            double eff_resist = dis[edge_point1] + dis[edge_point2] - 2*dis[LCA_point];
            edge.eff_w = eff_resist * edge.w;
            copy_off_tree_edge[i] = edge;
        }
    }
    printTime("resistance: E-V+1(off-tree) get_LCA")
}

void write_edge(vector<edge_t> &edge, const char *file){
    FILE *fp;
    fp = fopen(file, "w");
    for (int i=0; i<edge.size(); i++) {
        fprintf(fp, "%2d %2d %.20f \t%.20f\n", edge[i].u, edge[i].v, edge[i].eff_w, edge[i].w);
        // fout<<edge_point1<<" "<<edge_point2<<" "<<edge[i][2]<<" "<<edge[i][3]<<endl;
    }
    // fout.close();
    fclose(fp);
}

void write_bfs_process(vector<int> &bfs_process1, const char *file){
    FILE *fp;
    fp = fopen(file, "w");
    for (int i=0; i<bfs_process1.size(); i++) {
        fprintf(fp, "%d\n", bfs_process1[i]);
    }
    fclose(fp);
}
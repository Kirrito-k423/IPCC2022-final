#include "global.h"

vector<vector<vector<double>>> adja_list;
double *dis;
int *parent;
int *no_weight_dis;

/**
 * 根据树的边集表示，创建邻接表
*/
void create_adja_list(vector<vector<double>> &tree, vector<vector<vector<double>>> &adja_list){
    // vector<vector<int>> adja_list(N);
    int edge_num = tree.size();
    for(int i=0; i<edge_num; i++){
        int p1 = tree[i][0]-1;
        int p2 = tree[i][1]-1;
        // printf("%d %d %f\n", p1, p2, tree[i][3]);
        vector<double> elem;    //[point_idx, edge_weight]
        elem.push_back(p2);
        elem.push_back(tree[i][3]);
        adja_list[p1].push_back(elem);
        vector<double> elem2;    //[point_idx, edge_weight]
        elem2.push_back(p1);
        elem2.push_back(tree[i][3]);
        adja_list[p2].push_back(elem2);
    }
}

/**
 * DFS遍历生成树获得：
 * 1. 每个点到根节点的距离
 * 2. 每个点的parent数组
 * 3. 每个点到根节点的无权重距离
 * 根节点使用largest_volume_point
*/
void DFS_traversal(vector<vector<vector<double>>> &adja_list, double *dis, int *parent, int *no_weight_dis){
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
        for(int i=0; i<adja_list[point].size(); i++){
            int search_point = adja_list[point][i][0];
            if(no_weight_dis[search_point]==-1){   //first search
                process.push(search_point);
                dis[search_point] = dis[point] + 1.0/adja_list[point][i][1];    //add edge weight
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

void debug_print_adja(vector<vector<vector<double>>> &adja_list){
    printf("adja: \n");
    for(int i=0; i<M; i++){
        printf("%2d: ", i);
        for(int j=0; j<adja_list[i].size(); j++){
            printf("%d ", int(adja_list[i][j][0]));
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
void caculate_resistance(vector<vector<double>> &spanning_tree, vector<vector<double>> &off_tree_edge, vector<vector<double>> &copy_off_tree_edge){
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);

    adja_list.resize(M);
    dis = (double *)malloc(M*sizeof(double));   //每个点到根节点距离
    parent = (int *)malloc(M*sizeof(int));         //每个点父节点数组
    no_weight_dis = (int *)malloc(M*sizeof(int));  //每个点到根节点无权重距离

    //创建生成树的邻接表
    create_adja_list(spanning_tree, adja_list);
    printTime("adja list\t\t took %f ms\n")
    // debug_print_adja(adja_list);
    
    //DFS遍历，得到dis,parent等
    DFS_traversal(adja_list, dis, parent, no_weight_dis);
    printTime("DFS traversal\t\t took %f ms\n")
    // debug_print(dis, parent, no_weight_dis);
    
    // printf("path: \n");
    // debug_print_path(1859, 3044, parent, no_weight_dis);

    vector<double> edge;
    for (int i=0; i<off_tree_edge.size(); i++) {
        int edge_point1 = int(off_tree_edge[i][0])-1;
        int edge_point2 = int(off_tree_edge[i][1])-1;
        edge.push_back(edge_point1+1);
        edge.push_back(edge_point2+1);

        //树上计算等效电阻简单方式
        int LCA_point = get_LCA(edge_point1, edge_point2, parent, no_weight_dis);
        // printf("LCA(%d, %d)=%d\n", edge_point1, edge_point2, LCA_point);
        // int d_ = no_weight_dis[edge_point1] + no_weight_dis[edge_point2] - 2*no_weight_dis[LCA_point];
        // printf("no_weight_dis(%d, %d)=%d\n", edge_point1, edge_point2, d_);
        double eff_resist = dis[edge_point1] + dis[edge_point2] - 2*dis[LCA_point];
        edge.push_back(eff_resist * off_tree_edge[i][3]);

        edge.push_back(off_tree_edge[i][3]);
        copy_off_tree_edge.push_back(edge);
        edge.erase(edge.begin(),edge.end());
    }
    printTime("E-V+1(off-tree) get_LCA\t\t took %f ms\n")
}

void write_edge(vector<vector<double>> &edge, const char *file){
    FILE *fp;
    fp = fopen(file, "w");
    for (int i=0; i<edge.size(); i++) {
        int edge_point1 = int(edge[i][0]);
        int edge_point2 = int(edge[i][1]);
        fprintf(fp, "%2d %2d %.20f \t%.20f\n", edge_point1, edge_point2, edge[i][2], edge[i][3]);
        // fout<<edge_point1<<" "<<edge_point2<<" "<<edge[i][2]<<" "<<edge[i][3]<<endl;
    }
    // fout.close();
    fclose(fp);
}
#include "global.h"

bool compare(const edge_t &a,const edge_t &b){
    return a.eff_w > b.eff_w;
}

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


void kruscal(vector<edge_t> &edge_matrix, vector<edge_t> &spanning_tree){
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);                                                                         \
    int M = edge_matrix.size();
    //sort according to the weight of each edge
    // stable_sort(edge_matrix.begin(), edge_matrix.end(), compare);
    // __gnu_parallel::stable_sort(edge_matrix.begin(), edge_matrix.end(), compare);
    p_mergesort<edge_t>(edge_matrix, SORT_NUM_THREADS, cmp);
    printTime("kruscal: Sort G edge")

    // //run kruscal to get largest-effect-weight spanning tree
    // // MEWST = maximum-effective-weight spanning tree
    // int assistance_size=M; //int(edge_matrix.size());
    // int assistance[assistance_size+1];//check wether some points construct the circle
    // for (int i=0; i<=assistance_size; i++) {
    //     assistance[i]=i;
    // }

    // int k=0;//show how many trees have been add into
    // int tmin;
    // int tmax;
    // //kruscal
    // for (int i=0; i<edge_matrix.size(); i++) {
    //     int edge_point1 = int(edge_matrix[i][0]);
    //     int edge_point2 = int(edge_matrix[i][1]);
    //     if (assistance[edge_point1]!=assistance[edge_point2]){
    //         k++;
    //         spanning_tree.push_back(edge_matrix[i]);
    //         tmin = assistance[edge_point1]>=assistance[edge_point2] ?assistance[edge_point2]:assistance[edge_point1];
    //         tmax = assistance[edge_point1]< assistance[edge_point2] ?assistance[edge_point2]:assistance[edge_point1];
    //         for (int j=1; j<=assistance_size; j++) {
    //             if (assistance[j]==tmin){
    //                 assistance[j]=tmax;
    //             }
    //         }
    //     }
    //     if (k==M-1){
    //         break;
    //     }
    // }

    DSU dsu(M+1);
    int edge_cnt = 0;   //记录添加的边数，达到M-1时停止
    for (int i=0; i<edge_matrix.size(); i++) {
        int u = int(edge_matrix[i].u);
        int v = int(edge_matrix[i].v);
        if(!dsu.same_set(u, v)){
            spanning_tree.push_back(edge_matrix[i]);
            dsu.unite(u, v);
            edge_cnt++;
        }
        if(edge_cnt == M-1){
            break;
        }
    }
    printTime("kruscal: loop");
}

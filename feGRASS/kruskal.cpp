#include "global.h"

bool compare(const vector<double> &a,const vector<double> &b){
    return a[2]>b[2];
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


void kruscal(vector<vector<double>> &edge_matrix, vector<vector<double>> &spanning_tree){
    gettimeofday(&startTime, NULL);                                                                         \
    int M = edge_matrix.size();
    //sort according to the weight of each edge
    stable_sort(edge_matrix.begin(), edge_matrix.end(), compare);
    printTime("Sort G edge\t took %f ms\n")

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

    DSU dsu(M);
    int edge_cnt = 0;   //记录添加的边数，达到M-1时停止
    for (int i=0; i<edge_matrix.size(); i++) {
        int u = int(edge_matrix[i][0]);
        int v = int(edge_matrix[i][1]);
        if(!dsu.same_set(u, v)){
            spanning_tree.push_back(edge_matrix[i]);
            dsu.unite(u, v);
            edge_cnt++;
        }
        if(edge_cnt == M-1){
            break;
        }
    }
}

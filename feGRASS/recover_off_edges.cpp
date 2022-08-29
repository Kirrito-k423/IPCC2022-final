/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-29 19:59:51
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-08-29 21:17:11
 */
#include "global.h"


int calculate_belta(int i, MatrixXd *LG, int largest_volume_point, int edge_point1, int edge_point2){
    //run tarjan algorithm to get the upper bound
    stack<int> process;//to show the process of dfs
    int mark[M+1];//to show whether a point has been gone through
    int find[M+1];//Joint search set
    int position_node[M+1];//restore the child point in the dfs at the time (代表已经遍历到的位置，便于下一次恢复)
    for (int j=0; j<M+1; j++) {
        mark[j]=0;
        find[j]=j;
        position_node[j]=1;
    }
    process.push(largest_volume_point);//choose largest-volume point as root point
    mark[largest_volume_point]=1;

    //stop the search when the vertexes of edge have been found
    while (mark[edge_point1]==0||mark[edge_point2]==0) { //如果边的两点都没有遍历到
        for (int k=position_node[process.top()]-1; k<N; k++) { //第一次while。 position_node初始值全1，k初始值是0
            if ((* LG)((process.top())-1,k) != 0 &&    //k 和 top之间有边
                    mark[k+1] != 1 &&              //k 点没有遍历过
                    k != (process.top())-1 ) {      //k 和 top 不是同一个点
                position_node[process.top()]=k+1;   
                mark[k+1]=1;                        
                process.push(k+1);                  
                break;                              //mark 了一个点，就while check一下
            }
            else if (k==N-1){
                int a=process.top();
                process.pop();
                find[a]=process.top();              //貌似想构造以largest_volume_point为根的树，但是找到 两点the vertexes of edge又中断了
            }
        }
    }
    //get the first point we found in dfs
    int node=process.top();
    if (node==edge_point1) {
        node=edge_point2;
    }
    else {
        node=edge_point1;
    }

    DEBUG_PRINT("copy_off_tree_edge Loop %d/%ld \t first_point\t %d \t second_point\t %d\n",i,copy_off_tree_edge.size(),node,process.top());

    //attain the no-weight distance between the first point that we found and LCA
    int d1=0;
    while (true) {
        node=find[node];
        d1++;
        if (find[node]==node) {
            break;
        }
    }

    DEBUG_PRINT("copy_off_tree_edge Loop %d/%ld \t LCA node \t %d \t large_vol_point\t %d\n",i,copy_off_tree_edge.size(),node,largest_volume_point);

    //attain the no-weight distance between the second point and LCA
    int d2=1;
    while (true) {
        process.pop();
        if (process.top()==node) {
            break;
        }
        d2++;
    }

    DEBUG_PRINT("copy_off_tree_edge Loop %d/%ld \t d1 \t\t %d \t d2 \t\t %d\n",i,copy_off_tree_edge.size(),d1,d2);

    //compare between two nodes to get the lower one to limit the upper bound of bfs
    return d1>=d2?d2:d1;
}

void belta_BFS(int belta, MatrixXd *LG, std::vector<int> *candidate_point_set, int edge_point){
    (* candidate_point_set).push_back(edge_point);
    //use zero to cut the near layer
    (* candidate_point_set).push_back(0);
    int mark1[M+1];
    for (int j=0; j<M+1; j++) {
        mark1[j]=0;
    }
    mark1[edge_point]=1;
    int laywer=0;
    for (int j=0;j<M;j++) {
        if (laywer==belta){
            break;
        }
        if ((* candidate_point_set)[j]==0){
            (* candidate_point_set).push_back(0);
            laywer++;
        }
        else{
            for (int x=0; x<M; x++) {
                if ((* LG)((* candidate_point_set)[j]-1,x)==0) { //第一次是寻找和root节点相联的边
                    continue;
                }
                else if(mark1[x+1] == 0 &&          //没遍历过
                        x+1 != (* candidate_point_set)[j] &&   //不是自己
                        x+1 != edge_point ){       //不是root节点
                    (* candidate_point_set).push_back(x+1);
                    mark1[x+1]=1;
                }
            }
        }
    }
}

void adjust_similarity_tree(int i, std::vector<int> *bfs_process1, std::vector<int> *bfs_process2 ,\
                            int *similarity_tree, vector<vector<double>> *copy_off_tree_edge){
    //mark the edge that is similar to the edge which wants to be added
    for (int j=0; j<(* bfs_process1).size(); j++) {
        if ((* bfs_process1)[j]==0) {
            continue;
        }
        for (int k=0; k<(* bfs_process2).size(); k++) {
            if ((* bfs_process2)[k]==0) {
                continue;
            }
            if ((* bfs_process1)[j]==(* bfs_process2)[k]) {
                continue;
            }
            for (int z=i; z<(* copy_off_tree_edge).size(); z++) { // 余下的off_edge里，如果该边的两点，有一点在两个bfs的点集里，则该边视作similar
                if ((* copy_off_tree_edge)[z][0]==(* bfs_process1)[j]&&
                    (* copy_off_tree_edge)[z][1]==(* bfs_process2)[k]) {
                    similarity_tree[z]=1;
                }
                else if ((* copy_off_tree_edge)[z][0]==(* bfs_process2)[k]&&
                        (* copy_off_tree_edge)[z][1]==(* bfs_process1)[j]){
                    similarity_tree[z]=1;
                }
            }
        }
    }
}
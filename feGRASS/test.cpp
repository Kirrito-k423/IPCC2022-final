/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-31 20:16:21
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-08-31 22:52:15
 */


#include "global.h"

void test_LCA_find_update(MatrixXd *LG, int largest_volume_point){
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

    int marked_num=1;

    //stop the search when the vertexes of edge have been found
    while (marked_num < M) { //所有的点都要mark到
        for (int k=position_node[process.top()]-1; k<N; k++) { //第一次while。 position_node初始值全1，k初始值是0
            if ((* LG)((process.top())-1,k) != 0 &&    //k 和 top之间有边
                    mark[k+1] != 1 &&              //k 点没有遍历过
                    k != (process.top())-1 ) {      //k 和 top 不是同一个点
                position_node[process.top()]=k+1;   
                mark[k+1]=1;       
                marked_num++;                 
                process.push(k+1);                  
                break;                              //mark 了一个点，就while check一下
            }
            else if (k==N-1){
                int a=process.top();
                process.pop();
                if(find[a]==a){
                    DEBUG_PRINT("find[%d\t] update from %d \t to %d\n", a, find[a], process.top());
                }else{
                    DEBUG_PRINT("\t overwrite find[%d\t] update from %d \t to %d\n", a, find[a], process.top());
                }
                find[a]=process.top();              //貌似想构造以largest_volume_point为根的树，但是找到 两点the vertexes of edge又中断了
            }
        }
    }
    DEBUG_PRINT("pop ----------\n");
    while(!process.empty()){
        int a=process.top();
        process.pop();
        if(process.empty())
            break;
        if(find[a]==a){
            DEBUG_PRINT("find[%d\t] update from %d \t to %d\n", a, find[a], process.top());
        }else{
            DEBUG_PRINT("\t overwrite find[%d\t] update from %d \t to %d\n", a, find[a], process.top());
        }
        find[a]=process.top();              //貌似想构造以largest_volume_point为根的树，但是找到 两点the vertexes of edge又中断了
    }
}
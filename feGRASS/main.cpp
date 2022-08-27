#include "Eigen/Dense"
#include "Eigen/LU"
#include <stack>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <sys/time.h>
#include "global.h"

#define ROW 0
#define COLUMN 1
#define VALUE 2
#define vertex1 0
#define vertex2 1
#define Weight 2

using namespace Eigen;
using namespace std;

bool compare(const vector<double> &a,const vector<double> &b){
    return a[2]>b[2];
}

int main(int argc, const char * argv[]) {
    double startTime = omp_get_wtime();
    //read input file
    const char* file = "byn1.mtx";
    if(argc > 2) {
        printf("Usage : ./main <filename>");
        exit(-1);
    } else if(argc == 2) {
        file = argv[1];
    }
    //the matrix you read must be a adjacency matrix
    ifstream fin(file);
    int M, N, L;
    //Ignore headers and comments
    while (fin.peek() == '%')
        fin.ignore(2048, '\n');
    //declare matrix vector, volume and degree
    fin>>M>>N>>L;
    double volume[M+1];//volume of every point
    double degree[M+1];//degree of every point
    degree[0]=0;
    volume[0]=0;
    vector<double> triple;//element
    vector<vector<double>> triple1;//elements in the same column
    vector<vector<vector<double>>> triple2;//the whole matrix

    //fill matrix vector and calculate the degree and volume of every point
    int column=1;
    int m, n=0;
    double data=0;
    int n1=1;
    for (int i = 0; i < L; ++i) {
        fin>>m>>n>>data;
        //check wether change the column
        if (n1!=n){
            column++;
            //initialization
            volume[column]=0;
            degree[column]=0;
            n1=n;
            triple2.push_back(triple1);
            //clear up the triple1
            triple1.erase(triple1.begin(),triple1.end());
        }
        volume[column]=volume[column]+data;
        degree[column]=degree[column]+1;
        triple.push_back(m);
        triple.push_back(n);
        triple.push_back(data);
        triple1.push_back(triple);
        triple.erase(triple.begin(),triple.end());
    }
    fin.close();
    triple2.push_back(triple1);
    //free the memory of triple1
    vector<vector<double>>().swap(triple1);

    TIME_PRINT("Before timing \t took %f ms\n", 1000*(omp_get_wtime() - startTime));

    /**************************************************/
    /***************** Start timing *******************/
    /**************************************************/
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //find the point that has the largest volume
    int largest_volume_point=0;
    double largest_volume=0;
    for (int i=1;i<=M;i++){
        if (volume[i]>largest_volume){
            largest_volume_point=i;
            largest_volume=volume[i];
        }
    }

    //run bfs to get the no-weight distance between normal point with the largest-volume point
    int no_weight_distance[M+1];//no-weight-distance between largest-volume point and normal point
    for (int i=0;i<M+1;i++){
        no_weight_distance[i]=0;
    }
    queue<int> process;//to show the process of bfs
    //run bfs and calculate the no-weight distance
    int distance=1;
    int point;
    process.push(largest_volume_point);
    //use 0 to cut the layer
    process.push(0);
    no_weight_distance[largest_volume_point]=-1;
    while(process.size()!=1||process.front()!=0){
        point=process.front();
        process.pop();
        //set zero to cut near layers
        if (point==0){
            process.push(0);
            distance++;
            continue;
        }
        for(int i=0;i<triple2[point-1].size();i++){
            if(no_weight_distance[int(triple2[point-1][i][ROW])]==0){
                process.push(int(triple2[point-1][i][ROW]));
                no_weight_distance[int(triple2[point-1][i][ROW])]=distance;
            }
        }
    }
    no_weight_distance[largest_volume_point]=0;
    //free the memory
    queue<int>().swap(process);

    //construct the edge-weight matrix
    vector<double> edge;
    vector<vector<double>> edge_matrix;//to run the krascal
    for (int i=0; i<triple2.size(); i++) {
        for (int j=0; j<triple2[i].size(); j++) {
            if (triple2[i][j][COLUMN]<=triple2[i][j][ROW]){
                break;
            }
            edge.push_back(triple2[i][j][ROW]);
            edge.push_back(triple2[i][j][COLUMN]);
            //push effect valuable into vector
            double a=degree[int(edge[0])];
            double b=degree[int(edge[1])];
            double c=a>=b?a:b;
            double d=triple2[i][j][VALUE];
            double e=no_weight_distance[int(edge[0])];
            double f=no_weight_distance[int(edge[1])];
            double g=d*log(c)/(f+e);
            edge.push_back(g);
            edge.push_back(d);
            edge_matrix.push_back(edge);
            edge.erase(edge.begin(),edge.end());
        }
    }

    //sort according to the weight of each edge
    stable_sort(edge_matrix.begin(), edge_matrix.end(), compare);

    //run kruscal to get largest-effect-weight spanning tree
    int assistance[int(edge_matrix.size())+1];//check wether some points construct the circle
    for (int i=0; i<=edge_matrix.size(); i++) {
        assistance[i]=i;
    }
    int k=0;//show how many trees have been add into
    vector<vector<double>> spanning_tree;//spanning tree
    int tmin;
    int tmax;
    //kruscal
    for (int i=0; i<edge_matrix.size(); i++) {
        if (assistance[int(edge_matrix[i][0])]!=assistance[int(edge_matrix[i][1])]){
            k++;
            spanning_tree.push_back(edge_matrix[i]);
            tmin=assistance[int(edge_matrix[i][0])]>=assistance[int(edge_matrix[i][1])]?assistance[int(edge_matrix[i][1])]:assistance[int(edge_matrix[i][0])];
            tmax=assistance[int(edge_matrix[i][0])]<assistance[int(edge_matrix[i][1])]?assistance[int(edge_matrix[i][1])]:assistance[int(edge_matrix[i][0])];
            for (int j=1; j<=int(edge_matrix.size()); j++) {
                if (assistance[j]==tmin){
                    assistance[j]=tmax;
                }
            }
        }
        if (k==M-1){
            break;
        }
    }

    //construct the off-tree edge
    vector<vector<double>> off_tree_edge;
    int inside=0;//To show which trees are in the spanning tree
    for (int i=0; i<edge_matrix.size(); i++) {
        if (edge_matrix[i][0]==spanning_tree[inside][0]&&edge_matrix[i][1]==spanning_tree[inside][1]) {
            inside++;
            //to avoid inside crossing the border of spanning_tree
            if (inside==spanning_tree.size()) {
                inside--;
            }
            continue;
        }
        off_tree_edge.push_back(edge_matrix[i]);
    }
    vector<vector<double>>().swap(edge_matrix);

    //attain the Laplace matrix of spanning tree
    //restore in SparseMatrix of Eigen
    //that can get inverse matrix by LU decompose

    //construct the Laplace matrix of spanning tree
    MatrixXd LG=Eigen::MatrixXd::Zero(M,N);
    for (int i=0; i<spanning_tree.size(); i++) {
        LG(int(spanning_tree[i][0])-1,int(spanning_tree[i][1])-1)=-spanning_tree[i][3];
        LG(int(spanning_tree[i][1])-1,int(spanning_tree[i][0])-1)=-spanning_tree[i][3];
        LG(int(spanning_tree[i][1])-1,int(spanning_tree[i][1])-1)+=spanning_tree[i][3];
        LG(int(spanning_tree[i][0])-1,int(spanning_tree[i][0])-1)+=spanning_tree[i][3];
    }
    MatrixXd pseudo_inverse_LG=(LG.transpose()*LG).inverse()*LG.transpose();

    //calculate the resistance of each off_tree edge
    vector<vector<double>> copy_off_tree_edge;//to resore the effect resistance
    edge.erase(edge.begin(),edge.end());
    for (int i=0; i<off_tree_edge.size(); i++) {
            edge.push_back(off_tree_edge[i][0]);
            edge.push_back(off_tree_edge[i][1]);
            double a=pseudo_inverse_LG(int(edge[0])-1,int(edge[0])-1);
            double b=pseudo_inverse_LG(int(edge[1])-1,int(edge[1])-1);
            double c=pseudo_inverse_LG(int(edge[1])-1,int(edge[0])-1);
            double d=pseudo_inverse_LG(int(edge[0])-1,int(edge[1])-1);
            edge.push_back(a+b-c-d);
            edge.push_back(off_tree_edge[i][3]);
            copy_off_tree_edge.push_back(edge);
            edge.erase(edge.begin(),edge.end());
        }

    //sort by effect resistance
    vector<vector<double>>().swap(off_tree_edge);
    stable_sort(copy_off_tree_edge.begin(), copy_off_tree_edge.end(), compare);

    for (int i=0; i<M; i++) { if(i>30){break;}
        for (int j=0; j<N; j++) { if(j>30){break;}
            cout<<LG(i,j)<<' ';
        }
        cout<<endl;
    }
    //add some edge into spanning tree
    int num_additive_tree=0;
    int similarity_tree[copy_off_tree_edge.size()];//check whether a edge is similar to the edge added before
    for (int i=0; i<copy_off_tree_edge.size(); i++) {
        similarity_tree[i]=0;
    }
    for (int i=0; i<copy_off_tree_edge.size(); i++) {
        //if there has enough off-tree edge added into spanning tree, the work has been finished
        if (num_additive_tree==max(int(copy_off_tree_edge.size()/25), 2)) {
            break;
        }
        //if adge is not the similar tree,you can add it into spanning tree
        if (similarity_tree[i]==0){
            num_additive_tree++;
            /**** Iteration Log. Yuo delete the printf call. ****/
            if ((num_additive_tree%64)==0) {
                printf("num_additive_tree : %d\n", num_additive_tree);
            }
            spanning_tree.push_back(copy_off_tree_edge[i]);
            //run tarjhan algorithm to get the upper bound
            stack<int> process;//to show the process of dfs
            int mark[M+1];//to show whether a point has been gone through
            int find[M+1];//Joint search set
            int position_node[M+1];//restore the child point in the dfs at the time
            for (int j=0; j<M+1; j++) {
                mark[j]=0;
                find[j]=j;
                position_node[j]=1;
            }
            process.push(largest_volume_point);//choose largest-volume point as root point
            mark[largest_volume_point]=1;

            //stop the search when the vertexes of edge have been found
            while (mark[int(copy_off_tree_edge[i][0])]==0||mark[int(copy_off_tree_edge[i][1])]==0) {
                for (int k=position_node[process.top()]-1; k<N; k++) {
                    if (LG((process.top())-1,k)!=0&&mark[k+1]!=1&&k!=(process.top())-1) {
                        position_node[process.top()]=k+1;
                        mark[k+1]=1;
                        process.push(k+1);
                        break;
                    }
                    else if (k==N-1){
                        int a=process.top();
                        process.pop();
                        find[a]=process.top();
                    }
                }
            }
            //get the first point we found in dfs
            int node=process.top();
            if (node==copy_off_tree_edge[i][0]) {
                node=copy_off_tree_edge[i][1];
            }
            else {
                node=copy_off_tree_edge[i][0];
            }

            //attain the no-weight distance between the first point that we found and LCA
            int d1=0;
            while (true) {
                node=find[node];
                d1++;
                if (find[node]==node) {
                    break;
                }
            }

            //attain the no-weight distance between the second point and LCA
            int d2=1;
            while (true) {
                process.pop();
                if (process.top()==node) {
                    break;
                }
                d2++;
            }

            //compare between two nodes to get the lower one to limit the upper bound of bfs
            int belta=d1>=d2?d2:d1;

            //choose two nodes as root node respectively to run belta bfs
            vector<int> bfs_process1;
            bfs_process1.push_back(copy_off_tree_edge[i][0]);
            //use zero to cut the near layer
            bfs_process1.push_back(0);
            int mark1[M+1];
            for (int j=0; j<M+1; j++) {
                mark1[j]=0;
            }
            mark1[int(copy_off_tree_edge[i][0])]=1;
            int laywer=0;
            for (int j=0;j<M;j++) {
                if (laywer==belta){
                    break;
                }
                if (bfs_process1[j]==0){
                    bfs_process1.push_back(0);
                    laywer++;
                }
                else{
                    for (int x=0; x<M; x++) {
                        if (LG(bfs_process1[j]-1,x)==0) {
                            continue;
                        }
                        else if(mark1[x+1]==0&&x+1!=bfs_process1[j]&&x+1!=copy_off_tree_edge[i][0]){
                            bfs_process1.push_back(x+1);
                            mark1[x+1]=1;
                        }
                    }
                }
            }

            vector<int> bfs_process2;
            bfs_process2.push_back(copy_off_tree_edge[i][1]);
            //use zero to cut the near layer
            bfs_process2.push_back(0);
            int mark2[M+1];
            for (int j=0; j<M+1; j++) {
                mark2[j]=0;
            }
            mark2[int(copy_off_tree_edge[i][1])]=1;
            laywer=0;
            for (int j=0;j<M;j++) {
                if (laywer==belta){
                    break;
                }
                if (bfs_process2[j]==0){
                    bfs_process2.push_back(0);
                    laywer++;
                }
                else{
                    for (int x=0; x<M; x++) {
                        if (LG(bfs_process2[j]-1,x)==0) {
                            continue;
                        }
                        else if(mark2[x+1]==0&&x+1!=bfs_process2[j]&&x+1!=copy_off_tree_edge[i][1]){
                            bfs_process2.push_back(x+1);
                            mark2[x+1]=1;
                        }
                    }
                }
            }

            //mark the edge that is similar to the edge which wants to be added
            for (int j=0; j<bfs_process1.size(); j++) {
                if (bfs_process1[j]==0) {
                    continue;
                }
                for (int k=0; k<bfs_process2.size(); k++) {
                    if (bfs_process2[k]==0) {
                        continue;
                    }
                    if (bfs_process1[j]==bfs_process2[k]) {
                        continue;
                    }
                    for (int z=i; z<copy_off_tree_edge.size(); z++) {
                        if (copy_off_tree_edge[z][0]==bfs_process1[j]&&copy_off_tree_edge[z][1]==bfs_process2[k]) {
                            similarity_tree[z]=1;
                        }
                        else if (copy_off_tree_edge[z][0]==bfs_process2[k]&&copy_off_tree_edge[z][1]==bfs_process1[j]){
                            similarity_tree[z]=1;
                        }
                    }
                }
            }
        }
    }

    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
    /**************************************************/
    /******************* End timing *******************/
    /**************************************************/

    FILE* out = fopen("result.txt", "w");
    for(int i=0; i<spanning_tree.size(); i++){
        fprintf(out, "%d %d\n", int(spanning_tree[i][0]), int(spanning_tree[i][1]));
    }
    fclose(out);
    return 0;
}

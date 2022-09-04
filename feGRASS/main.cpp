#include "global.h"


#define ROW 0
#define COLUMN 1
#define VALUE 2
#define vertex1 0
#define vertex2 1
#define Weight 2

int M;
int N;
int L;
int largest_volume_point;
double subTime[5]={0,0,0,0,0}; // 伪逆， 循环总时间， 循环内三部分时间

bool compare(const vector<double> &a,const vector<double> &b){
    return a[2]>b[2];
}

void print_time_proportion(double total_time){
    printf("等效电阻\t循环总时间\t belta\t\t 2 BFS\t\t OMP_similarity\n");
    int i;
    int length=sizeof(subTime)/sizeof(subTime[0]);
    for(i=0; i<length-1; i++){
        printf("%8.2f\t",subTime[i]);
    }
    printf("%8.2f\n",subTime[i]);
    for(i=0; i<length-1; i++){
        printf("%4.2f%%\t\t",100*subTime[i]/total_time);
    }
    printf("%4.2f%%\n",100*subTime[i]/total_time);
    printf("循环+伪逆 占比 %.2f%%\n", 100*(subTime[0]+subTime[1])/total_time);
}

int main(int argc, const char * argv[]) {
    struct timeval startTime, endTime;
    gettimeofday(&startTime, NULL);

    //read input file
    const char* file = "byn1.mtx";
    if(argc > 2) {
        printf("Usage : ./main <filename>");
        exit(0);
    } else if(argc == 2) {
        file = argv[1];
    }
    //the matrix you read must be a adjacency matrix
    ifstream fin(file);
    // int M, N, L;
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

    printTime("Before timing \t\t\t took %f ms\n");

    /**************************************************/
    /***************** Start timing *******************/
    /**************************************************/
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //find the point that has the largest volume
    largest_volume_point=1;
    double largest_volume=volume[1];
    for (int i=2;i<=M;i++){
        if (volume[i]>largest_volume){
            largest_volume_point=i;
            largest_volume=volume[i];
        }
    }

    printTime("find largestPoint \t\t took %f ms\n");

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
            int search_point=int(triple2[point-1][i][ROW]);
            if(no_weight_distance[search_point]==0){
                process.push(search_point);
                no_weight_distance[search_point]=distance;
            }
        }
    }
    no_weight_distance[largest_volume_point]=0;
    //free the memory
    queue<int>().swap(process);

    printTime("BFS point-distance \t\t took %f ms\n")

    //construct the edge-weight matrix
    vector<double> edge; // [point_index1, point_index2, W_eff, W_ij]
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
            edge_matrix.push_back(edge); // size = L/2, 总边数
            edge.erase(edge.begin(),edge.end());
        }
    }

    //sort according to the weight of each edge
    stable_sort(edge_matrix.begin(), edge_matrix.end(), compare);

    printTime("Const Sort edge-weight matrix\t took %f ms\n")

    //run kruscal to get largest-effect-weight spanning tree
    // MEWST = maximum-effective-weight spanning tree
    int assistance_size=M; //int(edge_matrix.size());
    int assistance[assistance_size+1];//check wether some points construct the circle
    for (int i=0; i<=assistance_size; i++) {
        assistance[i]=i;
    }
    int k=0;//show how many trees have been add into
    vector<vector<double>> spanning_tree;//spanning tree
    int tmin;
    int tmax;
    //kruscal
    for (int i=0; i<edge_matrix.size(); i++) {
        int edge_point1 = int(edge_matrix[i][0]);
        int edge_point2 = int(edge_matrix[i][1]);
        if (assistance[edge_point1]!=assistance[edge_point2]){
            k++;
            spanning_tree.push_back(edge_matrix[i]);
            tmin = assistance[edge_point1]>=assistance[edge_point2] ?assistance[edge_point2]:assistance[edge_point1];
            tmax = assistance[edge_point1]< assistance[edge_point2] ?assistance[edge_point2]:assistance[edge_point1];
            for (int j=1; j<=assistance_size; j++) {
                if (assistance[j]==tmin){
                    assistance[j]=tmax;
                }
            }
        }
        if (k==M-1){
            break;
        }
    }

    printTime("Run kruscal\t\t\t took %f ms\n")

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

    printTime("Construct off-tree edge\t\t took %f ms\n")

    //calculate the resistance of each off_tree edge
    vector<vector<double>> copy_off_tree_edge;//to resore the effect resistance
    caculate_resistance(spanning_tree, off_tree_edge, copy_off_tree_edge);
    // write_edge(spanning_tree, "edge-spanning_tree.log");
    // write_edge(copy_off_tree_edge, "edge-copy_off_tree_edge.log");

    struct timeval Matrix_end_time;
    gettimeofday(&Matrix_end_time, NULL);
    subTime[0]=(Matrix_end_time.tv_sec-startTime.tv_sec)*1000+(Matrix_end_time.tv_usec-startTime.tv_usec)/1000.0;
    printTime("Calculate resistance\t\t took %f ms\n")

    //sort by effect resistance
    vector<vector<double>>().swap(off_tree_edge);
    stable_sort(copy_off_tree_edge.begin(), copy_off_tree_edge.end(), compare);
    // write_edge(copy_off_tree_edge, "edge-copy_off_tree_edge-sort.log");

    printTime("Sort & Init add edge \t\t took %f ms\n");

    /*
    construct global off-tree edges hash map, for an off-tree edge 
    with point with index i and j, hash key could be i << 16 & j
    Attention: point index should supposed less than 2^16.
    */
    map<uint32_t, uint16_t> off_tree_edge_map;
    //ofstream ou;
    //string path = "off_tree_edges.txt";
    //ou.open(path.c_str());
    for(int i = 0; i < copy_off_tree_edge.size()/cut_similarity_range; i++){
        // store a pair of directed edges to reduce lookup time
        uint32_t key1 = (uint32_t(copy_off_tree_edge[i][0]) << 16) | uint32_t(copy_off_tree_edge[i][1]);
        uint32_t key2 = (uint32_t(copy_off_tree_edge[i][1]) << 16) | uint32_t(copy_off_tree_edge[i][0]);
        // DEBUG_PRINT("node 1: %x, node2: %x, key: %x, value: %d \n", uint32_t(copy_off_tree_edge[i][0]), uint32_t(copy_off_tree_edge[i][1]), key1, uint16_t(i));
        off_tree_edge_map[key1] = uint16_t(i);
        off_tree_edge_map[key2] = uint16_t(i);
        //off_tree_edge_map.insert(pair<uint32_t, uint16_t>(key, i));        
        // char s[100];
        // sprintf(s, "node 1: %x, node2: %x, key: %x, value: %d \n", uint32_t(copy_off_tree_edge[i][0]), uint32_t(copy_off_tree_edge[i][1]), key1, uint16_t(i));
        // ou << s;
    }
    //ou.close();
    printTime("Construct off-tree edge hash map\t\t took %f ms\n")

    //add some edge into spanning tree
    int num_additive_tree=0;
    int similarity_tree_length=copy_off_tree_edge.size()/cut_similarity_range;
    int similarity_tree[similarity_tree_length];//check whether a edge is similar to the edge added before
    for (int i=0; i<similarity_tree_length; i++) {
        similarity_tree[i]=0;
    }


    
    
    // // 动态大小的 任务池
    // int task_list[task_pool_size];
    // int similarity_tree_list[task_pool_size * similarity_tree_length];
    // 并行初始化？
    // memset(similarity_tree_list, 0, sizeof(int) * task_pool_size * similarity_tree_length);
    
    // int i=0;
    int max_num_additive_tree = max(int(copy_off_tree_edge.size()/25), 2);


    struct timeval loop_begin_time, loop_end_time;
    double tmp_past_time;
    gettimeofday(&loop_begin_time, NULL);
    printTime("before while \t\t took %f ms\n");

    for (int i=0; i<copy_off_tree_edge.size(); i++) {
        //if there has enough off-tree edge added into spanning tree, the work has been finished
        if (num_additive_tree == max_num_additive_tree) {
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

            int edge_point1 = int(copy_off_tree_edge[i][0]);
            int edge_point2 = int(copy_off_tree_edge[i][1]);
            // int belta = calculate_belta(i, &LG ,largest_volume_point, edge_point1, edge_point2 );
            // int belta = calculate_belta_from_find(i, find, edge_point1, edge_point2 );
            int beta = calculate_beta(edge_point1, edge_point2);
            DEBUG_PRINT("\nbeta %d",beta);

            gettimeofday(&endTime, NULL);
            tmp_past_time=(endTime.tv_sec-startTime.tv_sec)*1000+(endTime.tv_usec-startTime.tv_usec)/1000.0;
            subTime[2] += tmp_past_time;
            TIME_PRINT("\ncopy_off_tree_edge belta %d/%ld \t took %f ms\n",i,copy_off_tree_edge.size(), tmp_past_time);
            gettimeofday(&startTime, NULL);

            //choose two nodes as root node respectively to run belta bfs
            vector<int> bfs_process1;
            beta_BFS(beta, bfs_process1, edge_point1);
            // printf("%d(%d): ", edge_point1, belta);
            // for(int j=0; j<bfs_process1.size(); j++){
            //     printf("%d ", bfs_process1[j]);
            // }
            // printf("\n");

            vector<int> bfs_process2;
            beta_BFS(beta, bfs_process2, edge_point2);

            gettimeofday(&endTime, NULL);
            tmp_past_time=(endTime.tv_sec-startTime.tv_sec)*1000+(endTime.tv_usec-startTime.tv_usec)/1000.0;
            subTime[3] += tmp_past_time;
            TIME_PRINT("copy_off_tree_edge 2 BFS %d/%ld \t took %f ms\n",i,copy_off_tree_edge.size(), tmp_past_time);
            gettimeofday(&startTime, NULL);

            DEBUG_PRINT("copy_off_tree_edge Loop %d/ \t bfs_process1 \t %ld \t bfs_process2\t %ld \t 3X \t %ld\n",i,
                        bfs_process1.size(),bfs_process2.size(),bfs_process1.size()*bfs_process2.size()*(copy_off_tree_edge.size()-i));

            // original adjust similarity tree
            // adjust_similarity_tree(i, bfs_process1, bfs_process2, similarity_tree, copy_off_tree_edge);

            // using hash map to store off tree edges
            // DEBUG_PRINT("start to adjust similarity tree\n");
            adjust_similarity_tree(i, bfs_process1, bfs_process2, similarity_tree, off_tree_edge_map);

            // 假如按照论文，可以写同一个similarity_tree_list（不行，尝试过了，结果有几个是错的）
            gettimeofday(&endTime, NULL);
            tmp_past_time=(endTime.tv_sec-startTime.tv_sec)*1000+(endTime.tv_usec-startTime.tv_usec)/1000.0;
            subTime[4] += tmp_past_time;
            TIME_PRINT("copy_off_tree_edge similarity %d/%ld \t took %f ms\n",i,copy_off_tree_edge.size(), tmp_past_time);
            gettimeofday(&startTime, NULL);
        } 
    }

    gettimeofday(&loop_end_time, NULL);
    subTime[1]=(loop_end_time.tv_sec-loop_begin_time.tv_sec)*1000+(loop_end_time.tv_usec-loop_begin_time.tv_usec)/1000.0;
    // subTime[4] += (loop_end_time.tv_sec-startTime.tv_sec)*1000+(loop_end_time.tv_usec-startTime.tv_usec)/1000.0;
    printTime("\ncopy_off_tree_edge End \t\t took %f ms\n\n")

    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0);
    double total_time=(end.tv_sec-start.tv_sec)*1000+(end.tv_usec-start.tv_usec)/1000.0;
    /**************************************************/
    /******************* End timing *******************/
    /**************************************************/
    print_time_proportion(total_time);

    FILE* out = fopen("result.txt", "w");
    for(int i=0; i<spanning_tree.size(); i++){
        fprintf(out, "%d %d\n", int(spanning_tree[i][0]), int(spanning_tree[i][1]));
    }
    fclose(out);
    return 0;
}

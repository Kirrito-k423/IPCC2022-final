
#include "global.h"

typedef enum { FALSE,
               TRUE } boolean;

void parse_input(vector<vector<double>> &edges, int rank, int number_of_edges, int vertices_per_process, int edges_per_process);
void parse_edge_list_input();
int compare_edges(const void *a, const void *b);
boolean should_send();
int get_merge_partner_rank(boolean should_send, int rank, int merge_iteration);
void merge(int rank, int number_of_vertices, vector<vector<double>> &local_msf_edges, vector<vector<double>> &recv_msf_edges, vector<vector<double>> &merged_msf_edges);
void send_recieve_local_msf(int rank, int count, vector<vector<double>> &local_msf_edges, vector<vector<double>> & recv_msf_edges);
void merge_msf(int num_of_processors, int rank, int &merge_iteration);
void print_received_msf_edges();

std::ifstream fin;


typedef struct node_t { /* Union-Find data structure */
    struct node_t *parent;
    int depth;
} u_node;

u_node *uf_set; /* Array indicating which set does vertex belong to
                   (used in Union-Find algorithm) */

void uf_make(int number_of_vertices);
u_node *uf_find(u_node *a);
void uf_union(u_node *a, u_node *b);


bool compare(const vector<double> &a,const vector<double> &b){
    return a[2]<b[2];
}
boolean should_send(int rank, int merge_iteration) {
    if ((rank / (int)pow(2, merge_iteration)) % 2 == 0) {
        return FALSE;
    } else {
        return TRUE;
    }
}
int get_merge_partner_rank(boolean should_send, int rank, int merge_iteration) {
    // Merge partner is rank +/- 2^merge_iteration
    if (should_send) {
        return rank - pow(2, merge_iteration); // rank of destination
    } else {
        return rank + pow(2, merge_iteration); // rank of source
    }
}
void parse_input(vector<vector<double>> &edges, int rank, int number_of_edges, int vertices_per_process, int edges_per_process) {
    /* Start and end index of vertex belonging to this process */
    int rank_range_start = rank * vertices_per_process;
    int rank_range_end = (rank + 1) * vertices_per_process;
    printf("rank: %d\trange start:%d\trange end:%d\n", rank, rank_range_start, rank_range_end);

    int i;
    double f_edge[3];
    vector<double> edge;
    int u, v;
    double edge_weight;
    for (i = 0; i < number_of_edges; i++) {
        fin >> u >> v;
        fin >> edge_weight;
        //printf("%d\t%d\t%lf\n", (int)u, (int)v, edge_weight);
        if ((((int)u >= rank_range_start) && ((int)v < rank_range_end)) || (((int)u >= rank_range_start) && ((int)u < rank_range_end))) {
            edge.push_back((double)u);
            edge.push_back((double)v);
            edge.push_back(edge_weight);
            edges.push_back(edge);
            edge.erase(edge.begin(), edge.end());
        }
    }

    fin.close();
    MPI_Barrier(MPI_COMM_WORLD);
}


/*
每个进程根据分配的节点及与这些节点相连的边，构造局部生成树
local_msf_edge_count = merged_msf_edge_count = 局部生成树的边数
*/

void find_local_msf(int rank, int number_of_vertices, vector<vector<double>> & edges, vector<vector<double>> & local_msf_edges, vector<vector<double>> &merged_msf_edges) {

    /* Sort edges belonging to each process */
    stable_sort(edges.begin(), edges.end(), compare);
    //qsort(edges, number_of_edges, sizeof(edge_s), compare_edges);

    uf_make(number_of_vertices);
    vector<double> min_edge;
    int used_edge_index = 0;
    int i;
    for (i = 1; i <= edges.size(); i++) {
        min_edge = edges[used_edge_index++];

        int v = (int)min_edge[0];
        int u = (int)min_edge[1];
        u_node *v_node = uf_find(uf_set + v);
        u_node *u_node = uf_find(uf_set + u);

        /* Add edge to MSF if it doesn't form a cycle */
        if (v_node != u_node) {
            local_msf_edges.push_back(min_edge);
            merged_msf_edges.push_back(min_edge);
            uf_union(v_node, u_node);
        }
    }
}

/*
每个进程与其在迭代中的partner进程通信，
发送方：先发送边数，然后发送自己的局部生成树边集
接收方：接受边数，然后接收对方的局部生成树边集，然后执行merge
*/
void merge_msf(int number_of_processors, int number_of_vertices, int rank, vector<vector<double>> &local_msf_edges, vector<vector<double>> &merged_msf_edges) {
    int max_merge_iteration = 0;

    int temp = number_of_processors;
    while (temp > 1) {
        temp = temp / 2;
        max_merge_iteration++;
    }

    vector<vector<double>> recv_msf_edges;

    boolean has_sent = FALSE;
    int merge_iteration = 0;
    while (merge_iteration < max_merge_iteration) {
        MPI_Status status;
        printf("merge iteration: %d, max merge iteration: %d\n", merge_iteration, max_merge_iteration);
        boolean sending = should_send(rank, merge_iteration);
        int merge_partner_rank = get_merge_partner_rank(sending, rank, merge_iteration);
        if (sending) { // should send
            int local_msf_edges_size = local_msf_edges.size();
            // 先发送边数
            MPI_Send(&local_msf_edges_size, 1, MPI_INT, merge_partner_rank, 0, MPI_COMM_WORLD);
            MPI_Send(local_msf_edges.data(), local_msf_edges_size, MPI_DOUBLE, merge_partner_rank, 0, MPI_COMM_WORLD);
            printf("Proc %d: Sending %d edges to processor #%d\n", rank, local_msf_edges_size, merge_partner_rank);
            has_sent = TRUE;
        } else { // should receive
            int recv_msf_edges_size;
            MPI_Recv(&recv_msf_edges_size, 1, MPI_INT, merge_partner_rank, 0, MPI_COMM_WORLD, &status);
            recv_msf_edges.resize(recv_msf_edges_size);
            double *recv_data =  (double *) malloc(recv_msf_edges_size * sizeof(double) * 3);
            MPI_Recv(recv_data, recv_msf_edges_size, MPI_DOUBLE, merge_partner_rank, 0, MPI_COMM_WORLD, &status);
            vector<double> edge;
            for(int k = 0; k < recv_msf_edges_size; k++){
                edge.push_back(recv_data[k * 3]);
                edge.push_back(recv_data[k * 3 + 1]);
                edge.push_back(recv_data[k * 3] + 2);
                recv_msf_edges.push_back(edge);
                edge.erase(edge.begin(), edge.end());
            }
            //MPI_Get_count(&status, mpi_edge, &recv_msf_edge_count);
            merge(rank, number_of_vertices, local_msf_edges, recv_msf_edges, merged_msf_edges);
            printf("Proc %d: Received %d edges from processor #%d\n", rank, recv_msf_edges_size, merge_partner_rank);
        }
        merge_iteration++;
        if (has_sent) { // Process has sent it's data and can terminate now
            printf("Proc %d has send in iteration %d\n", rank, merge_iteration);
            break;
        }
        printf("rank: %d, iteration finished!\n", rank);
    }
    if(rank == 0){
        for (int i = 0; i < merged_msf_edges.size(); i++){
            printf("%d %d %lf\n", (int)merged_msf_edges[i][0], (int)merged_msf_edges[i][1], merged_msf_edges[i][2]);
        }
    }
}


/*
接收进程将接受到的局部生成树边，加入到自己的merged_msf_edges（上一次迭代产生的局部生成树边集）中，
然后将自己的merged_msf_edges边排序，从中不断选择权重较小的边，不产生回路时加入自己的局部生成树
*/
void merge(int rank, int number_of_vertices, vector<vector<double>> &local_msf_edges, vector<vector<double>> &recv_msf_edges, vector<vector<double>> &merged_msf_edges) {
    int i;
    printf("rank: %d, local size: %ld, merged size: %ld, recvd size: %ld\n", rank, local_msf_edges.size(), merged_msf_edges.size(), recv_msf_edges.size());
    for (i = 0; i < recv_msf_edges.size(); i++) {
        merged_msf_edges.push_back(recv_msf_edges[i]);
    }
    printf("rank: %d, local size: %ld, merged size: %ld, recvd size: %ld\n", rank, local_msf_edges.size(), merged_msf_edges.size(), recv_msf_edges.size());

    /* Sort local and received edges */
    //qsort(merged_msf_edges, merged_msf_edge_count, sizeof(), compare_edges);
    // for (int i = 0; i < merged_msf_edges.size(); i++){
    //     printf("%d %d %lf\n", (int)merged_msf_edges[i][0], (int)merged_msf_edges[i][1], merged_msf_edges[i][2]);
    // }
    stable_sort(merged_msf_edges.begin(), merged_msf_edges.end(), compare);
    uf_make(number_of_vertices);
    int used_edge_index = 0;
    vector<double> min_edge;
    for (i = 0; i < merged_msf_edges.size(); i++) {
        min_edge = merged_msf_edges[used_edge_index++];

        int v = (int)(min_edge[0]);
        int u = (int)(min_edge[1]);
        u_node *v_node = uf_find(uf_set + v);
        u_node *u_node = uf_find(uf_set + u);

        if (v_node != u_node) {
            local_msf_edges.push_back(min_edge);
            uf_union(v_node, u_node);
        }
    }

    /* Transfer new local MSF edges to merged edges */
    // merged_msf_edges.erase(merged_msf_edges.begin(), merged_msf_edges.end());
    // for (i = 0; i < local_msf_edges.size(); i++) {
    //     merged_msf_edges.push_back(local_msf_edges[i]);
    // }
    merged_msf_edges.erase(merged_msf_edges.begin(), merged_msf_edges.end());
    merged_msf_edges.assign(local_msf_edges.begin(), local_msf_edges.end());
}

void uf_make(int number_of_vertices) {
    int size = sizeof(u_node) * number_of_vertices; // + sizeof(int)*(number_of_vertices - 1);
    uf_set = (u_node *)malloc(size);
    memset(uf_set, 0, size);
}

u_node *uf_find(u_node *a) {
    if (a->parent == NULL)
        return a;
    else
        return (a->parent = uf_find(a->parent));
}

void uf_union(u_node *a, u_node *b) {
    if (a->depth > b->depth) {
        b->parent = a;
    } else if (a->depth < b->depth) {
        a->parent = b;
    } else {
        a->parent = b;
        a->depth += 1;
    }
}



int main(int argc, char **argv) {
    if (argc < 2) {
        printf("usage: %s filename\n", argv[0]);
        return 1;
    }
    printf("init mpi...\n");
    int rank, number_of_processors;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const char *file = "byn.mtx";
    if (argc > 2) {
        printf("Usage : ./main <filename>");
        exit(0);
    } else if (argc == 2) {
        file = argv[1];
    }
    int number_of_vertices, number_of_edges;
    fin = ifstream(file);
    while (fin.peek() == '%')
        fin.ignore(2048, '\n');
    fin >> number_of_vertices;
    fin >> number_of_vertices;
    fin >> number_of_edges;
    vector<vector<double>> edges; // 每个进程分配到的边
    vector<vector<double>> local_msf_edges;
    vector<vector<double>> merged_msf_edges;
    int vertices_per_process = number_of_vertices / number_of_processors;
    int edges_per_process = number_of_edges / number_of_processors;
    printf("parse input...\n");
    parse_input(edges, rank, number_of_edges, vertices_per_process, edges_per_process);

    printf("check conditions...\n");
    // check_conditions();

    printf("start parse measure...\n");
    // start_parse_measure();
    printf("end parse measure...\n");
    // end_parse_measure();

    printf("start computation measure...\n");
    printf("find local msf...\n");
    find_local_msf(rank, number_of_vertices, edges, local_msf_edges, merged_msf_edges);
    printf("merge msf...\n");
    merge_msf(number_of_processors, number_of_vertices, rank, local_msf_edges, merged_msf_edges);
    printf("end computation...\n");
    printf("print measurement times...\n");

    MPI_Finalize();

    return 0;
}
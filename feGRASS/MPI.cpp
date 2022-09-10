#include"global.h"

void print_mpi_vector(int *recv_buf, int length){
    for(int i=0; i<length; i++){
        DEBUG_PRINT("%d\t",recv_buf[i]);
        if(i%10==9)
            DEBUG_PRINT("\n");
    }
}

void fg_MPI_synchronization(vector<vector<int>> &syn_vector_list, int *similarity_tree){
    //统计各自要发送的数据个数
    int total_send_num=0;
    for (int k = 0; k < NUM_THREADS; k++) {
        total_send_num += syn_vector_list[k].size();
    }
    MPI_DEBUG_PRINT("mpi_rank %d\t total_send_num\t%d\n",mpi_rank,total_send_num);

    // allgather统计数据总数，声明空间
    int *vector_size_recvbuf = (int *)malloc(comm_size * sizeof(int));
    int *vector_size_displs = (int *)malloc(comm_size * sizeof(int));

    MPI_Allgather(&total_send_num, 1 ,MPI_INT,
                vector_size_recvbuf, 1, MPI_INT, MPI_COMM_WORLD);
    int All_vector_size=0;
    int displs=0;
    for(int i=0; i<comm_size; i++){
        All_vector_size += vector_size_recvbuf[i];
        vector_size_displs[i]=displs;
        displs += vector_size_recvbuf[i];
    }
    MPI_DEBUG_PRINT("mpi_rank %d\t All_vector_size\t%d\n",mpi_rank,All_vector_size);

    // allgatherv 接受数据
    int *send_buf = (int *)malloc(total_send_num * sizeof(int));
    int *recv_buf = (int *)malloc(All_vector_size * sizeof(int));
    for (int i = 0,k = 0; k < NUM_THREADS; k++) {
        for (int j = 0; j < syn_vector_list[k].size(); j++) {
            send_buf[i++] = syn_vector_list[k][j];
        }
    }
    MPI_Allgatherv(send_buf, total_send_num, MPI_INT,
                recv_buf, vector_size_recvbuf, vector_size_displs, MPI_INT, MPI_COMM_WORLD);
    // print_mpi_vector(recv_buf, All_vector_size);

    // 处理数据到similarity_tree里
    for (int k = 0; k < All_vector_size; k++) {
        similarity_tree[recv_buf[k]]=1;
    }

    free(vector_size_recvbuf);
    free(vector_size_displs);
    free(send_buf);
    free(recv_buf);
}

int* MPI_synchronization(int *vector_size_list, int *vector_displs_list,int vector_start, int vector_end, vector<vector<int>> &syn_vector_list){
    //统计各自要发送的数据个数
    int total_send_num = 0;
    int vector_num_sendbuf_size = (vector_end - vector_start);
    int *vector_num_sendbuf = (int *)malloc(vector_num_sendbuf_size * sizeof(int));
    for (int i = 0,k = vector_start; k < vector_end; k++) {
        vector_num_sendbuf[i++] = syn_vector_list[k].size();
        total_send_num += syn_vector_list[k].size();
    }
    MPI_DEBUG_PRINT("mpi_rank %d\t total_send_num\t%d\n",mpi_rank,total_send_num);

    // allgather统计数据总数，声明空间
    int *vector_size_recvbuf = (int *)malloc(comm_size* sizeof(int));
    int *vector_size_displs = (int *)malloc(comm_size * sizeof(int));

    MPI_Allgather(vector_num_sendbuf, vector_num_sendbuf_size ,MPI_INT,
                vector_size_list, vector_num_sendbuf_size, MPI_INT, MPI_COMM_WORLD);
    int All_vector_size=0;
    int displs=0;
    int index=0;
    for(int i=0; i<comm_size; i++){
        int tmp_rank_sum=0;
        vector_size_displs[i]=displs;
        for(int j=0;j<vector_num_sendbuf_size; j++){
            vector_displs_list[index]=displs;
            displs += vector_size_list[index];
            tmp_rank_sum += vector_size_list[index];
            index++;
            // DEBUG_PRINT("mpi_rank %d\t index tmp\t%d\t%d\n",mpi_rank,index,tmp_rank_sum);
        }
        // MPI_DEBUG_PRINT("mpi_rank %d\t tmp_rank_sum\t%d\n",mpi_rank,tmp_rank_sum);
        vector_size_recvbuf[i] = tmp_rank_sum;
        MPI_DEBUG_PRINT("mpi_rank %d\t buf displs\t%d\t%d\n",mpi_rank, vector_size_recvbuf[i],vector_size_displs[i]);
        All_vector_size += tmp_rank_sum; 
    }
    MPI_DEBUG_PRINT("mpi_rank %d\t All_vector_size\t%d\n",mpi_rank,All_vector_size);

    // allgatherv 接受数据
    int *send_buf = (int *)malloc(total_send_num * sizeof(int));
    int *recv_buf = (int *)malloc(All_vector_size * sizeof(int));
    for (int i = 0,k = vector_start; k < vector_end; k++) {
        for (int j = 0; j < syn_vector_list[k].size(); j++) {
            send_buf[i++] = syn_vector_list[k][j];
        }
    }
    MPI_Allgatherv(send_buf, total_send_num, MPI_INT,
                recv_buf, vector_size_recvbuf, vector_size_displs, MPI_INT, MPI_COMM_WORLD);
    // print_mpi_vector(recv_buf, All_vector_size);

    // 处理数据到syn_vector_list里

    free(vector_num_sendbuf);
    free(vector_size_recvbuf);
    free(vector_size_displs);
    free(send_buf);
    return recv_buf;
    // free(recv_buf);
}
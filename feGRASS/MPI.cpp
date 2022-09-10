#include"global.h"

void print_mpi_vector(int *recv_buf, int length){
    for(int i=0; i<length; i++){
        DEBUG_PRINT("%d\t",recv_buf[i]);
        if(i%10==9)
            DEBUG_PRINT("\n");
    }
}

void MPI_synchronization(vector<vector<int>> &syn_vector_list, int *similarity_tree){
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
    free(send_buf);
    free(recv_buf);
}
#include "global.h"

/**
 * caculate effect resistance of each off-tree edges on spanning tree
 * store result in copy_off_tree_edge
 * 
 * spanning_tree: [point_index1, point_index2, W_eff, W_ij]
 * off_tree_edge: subset of spanning_tree
 * copy_off_tree_edge: [point_index1, point_index2, R_ij, W_ij]
*/
void caculate_resistance(vector<vector<double>> &spanning_tree, vector<vector<double>> &off_tree_edge, vector<vector<double>> &copy_off_tree_edge, MatrixXd &LG){
    //MatrixXd pseudo_inverse_LG=(LG.transpose()*LG).inverse()*LG.transpose();
    MatrixXd pseudo_inverse_LG=LG.completeOrthogonalDecomposition().pseudoInverse();

    vector<double> edge; // [point_index1, point_index2, W_eff, W_ij]
    edge.erase(edge.begin(),edge.end());
    for (int i=0; i<off_tree_edge.size(); i++) {
        int edge_point1 = int(off_tree_edge[i][0]);
        int edge_point2 = int(off_tree_edge[i][1]);
        edge.push_back(edge_point1);
        edge.push_back(edge_point2);
        double a=pseudo_inverse_LG(edge_point1-1,edge_point1-1);
        double b=pseudo_inverse_LG(edge_point2-1,edge_point2-1);
        double c=pseudo_inverse_LG(edge_point2-1,edge_point1-1);
        double d=pseudo_inverse_LG(edge_point1-1,edge_point2-1);
        edge.push_back(a+b-c-d);
        edge.push_back(off_tree_edge[i][3]);
        copy_off_tree_edge.push_back(edge);
        edge.erase(edge.begin(),edge.end());
    }
}

void write_edge(vector<vector<double>> &edge, string file){
    ofstream fout(file);
    for (int i=0; i<edge.size(); i++) {
        int edge_point1 = int(edge[i][0]);
        int edge_point2 = int(edge[i][1]);
        fout<<edge_point1<<" "<<edge_point2<<" "<<edge[i][2]<<" "<<edge[i][3]<<endl;
    }
    fout.close();
}
#include "platform_atomics.hpp"
#include "pvector.h"
#include "sliding_queue.h"
#include "graph.h"
#include "timer.h"
#include "bitmap.h"
#include "global.h"


using namespace std;



void TDStep(pvector<NodeID> &parent,
               SlidingQueue<NodeID> &queue) {
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue);
    #pragma omp for nowait
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      NodeID u = *q_iter;
      for (edge_t v :adja_list[u]) {
        NodeID curr_val = parent[v.u];
        if (curr_val < 0) {
          if (compare_and_swap(parent[v.u], curr_val, u)) {
            lqueue.push_back(v.u);
          }
        }
      }
    }
    lqueue.flush();
  }
}

pvector<NodeID> InitParent() {
  pvector<NodeID> parent(M);
  #pragma omp parallel for
  for (NodeID n=0; n < M; n++)
    parent[n] = -1;
    // parent[n] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  return parent;
}


// g ==  vector<vector<edge_t>> adja_list;
// source = largest_volume_point 第一次
pvector<NodeID> DOBFS(NodeID source, int alpha = 15,
                      int beta = 18) {

  pvector<NodeID> parent = InitParent();
  parent[source] = source;
  SlidingQueue<NodeID> queue(M);
  queue.push_back(source);
  queue.slide_window();
  while (!queue.empty()) {
      TDStep(parent, queue);
      queue.slide_window();
  }
  #pragma omp parallel for
  for (NodeID n = 0; n < M; n++)
    if (parent[n] < -1)
      parent[n] = -1;
  return parent;
}
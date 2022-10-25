#ifndef __OUTPUT_GRAPH_H
#define __OUTPUT_GRAPH_H

#include <stdio.h>
#include <vector>
#include <map>
#include <Graph.h>

using namespace std;

void output_graph_yott(
    vector<tuple3> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    const int idx, 
    const int EDGE_SIZE,
    Graph &g
){
    
    int port[g.num_nodes()];
    for (int i = 0; i < g.num_nodes(); ++i) port[i] = 0;

    int a, b, dir, value;
    printf("digraph g {\n");

    for (int j = 0; j < g.num_nodes(); ++j) {
        printf("%d [label=\"%s\";op=\"%s\";]\n", j, g.get_name_node(j).c_str(), g.get_name_op(j).c_str());
    }

    for (int j = 0; j < EDGE_SIZE; ++j) {
        a = vector_edges[idx*EDGE_SIZE+j].v0;
        b = vector_edges[idx*EDGE_SIZE+j].v1;
        dir = vector_edges[idx*EDGE_SIZE+j].v2;
        if (dir == 0) printf("%d -> %d [port=%d", a, b, port[b]++);
        else printf("%d -> %d [port=%d", b, a, port[a]++);

        value = edges_cost[idx][make_pair(a,b)];
        if (value > 1) {
            printf(" color=red, label=%d", value);
        }
        printf("]\n");
    }
    printf("}\n");
} 

void output_graph_sa(
    vector<pair<int,int>> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    const int idx, 
    const int EDGE_SIZE, 
    Graph &g
){

    int port[g.num_nodes()];
    for (int i = 0; i < g.num_nodes(); ++i) port[i] = 0;
    
    int a, b, dir, value;
    printf("digraph g {\n");

    for (int j = 0; j < g.num_nodes(); ++j) {
        printf("%d [label=\"%s\";op=\"%s\";]\n", j, g.get_name_node(j).c_str(), g.get_name_op(j).c_str());
    }

    for (int j = 0; j < EDGE_SIZE; ++j) {
        a = vector_edges[j].first;
        b = vector_edges[j].second;
        
        printf("%d -> %d [port=%d", a, b, port[b]++);

        value = edges_cost[idx][make_pair(a,b)];
        if (value > 1) {
            printf(" color=red, label=%d", value);
        }
        printf("]\n");
    }
    printf("}\n");
}

#endif
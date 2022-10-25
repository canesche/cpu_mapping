#ifndef __STATISTIC_H
#define __STATISTIC_H

#include <Graph.h>

void statistic(Graph g, string name, int times, int *type_node, 
    vector<tuple3> vector_edges, vector<map<pair<int,int>,int>> edges_cost, 
    int *vector_cost) {

    const int NODE_SIZE = g.num_nodes();
    const int EDGE_SIZE = g.num_edges();
    const int INPUTS_SIZE = g.get_inputs().size();
    const int OUTPUTS_SIZE = g.get_outputs().size();

    int best_cost = 999999, index = 0;
    for (int i = 0; i < times; ++i) {
        if (best_cost > vector_cost[i]) {
            best_cost = vector_cost[i];
            index = i;
        }
    }

    float count_1 = 0, count_m_1 = 0;
    
    int total_edge = 0, worst_edge = -1, value;
    pair<int,int> key;
    for (int i = 0; i < EDGE_SIZE; ++i) {
        key = make_pair(vector_edges[i+index*EDGE_SIZE].v0, vector_edges[i+index*EDGE_SIZE].v1);
        value = edges_cost[index][key];
        total_edge += value;
        if (value > worst_edge) worst_edge = value;
        if (value == 1) count_1++;
        else count_m_1++;
    }
    count_1 = count_1*100 / EDGE_SIZE;
    count_m_1 = count_m_1*100 / EDGE_SIZE;

    int count_ff = 0, count_ioe = 0, count_rf = 0, count_w_ff = 0, count_w_rf = 0;

    int VISITED[NODE_SIZE];
    for (int i = 0; i < NODE_SIZE; ++i) VISITED[i] = 0;

    vector<pair<int,int>> e = g.get_edges();
    for (int i = 0; i < EDGE_SIZE; ++i) {
        int a = e[i].first;
        int b = e[i].second;
        
        if (VISITED[b] == 1) {
            count_w_rf += 1;
        } else {
            count_w_ff += 1;
            VISITED[a] = 1;
        }
    }
    
    for (int i = 0; i < NODE_SIZE; ++i) VISITED[i] = 0;

    e = g.get_edges();
    for (int i = 0; i < EDGE_SIZE; ++i) {
        int a = vector_edges[i+index*EDGE_SIZE].v0;
        int b = vector_edges[i+index*EDGE_SIZE].v1;
        
        if (type_node[a] == 1 || type_node[b] == 1) count_ioe += 1;
        else if (VISITED[b] == 1) {
            count_rf += 1;
        } else {
            count_ff += 1;
            VISITED[a] = 1;
        }
    }

    int count_ff_after = 0, count_ioe_after = 0, count_rf_after = 0;
    
    for (int i = 0; i < NODE_SIZE; ++i) VISITED[i] = 0;
    
    for (int i = 0; i < EDGE_SIZE; ++i) {
        int a = vector_edges[i+index*EDGE_SIZE].v0;
        int b = vector_edges[i+index*EDGE_SIZE].v1;
        
        key = make_pair(a,b);

        int new_cost = 1;
        if (edges_cost[index].count(key)) {
            new_cost = edges_cost[index][key];
        } else {
            key = make_pair(b, a);
            new_cost = edges_cost[index][key];
        }
        
        if (type_node[a] == 1 || type_node[b] == 1) count_ioe_after += new_cost;
        else if (VISITED[b] == 1) {
            count_rf_after += new_cost;
        } else {
            count_ff_after += new_cost;
            VISITED[a] = 1;
        }
    }

    printf("%s,%d,%d,%d/%d,%.2f,%d,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,", 
        name.c_str(), NODE_SIZE, EDGE_SIZE, INPUTS_SIZE, OUTPUTS_SIZE, 1.0*best_cost/EDGE_SIZE, worst_edge, 
        100.0*count_ff/EDGE_SIZE, 
        100.0*count_ioe/EDGE_SIZE, 
        100.0*count_rf/EDGE_SIZE, 
        100.0*count_w_ff/EDGE_SIZE, 
        100.0*count_w_rf/EDGE_SIZE, 
        1.0*count_ff_after/max(1,count_ff), 
        1.0*count_ioe_after/max(1,count_ioe),
        1.0*count_rf_after/max(1,count_rf),
        count_1,
        count_m_1
    ); 

}

#endif
#ifndef __TRY_ADJACENCY_H
#define __TRY_ADJACENCY_H

#include <vector>
#include <placed.h>
#include <get_cost.h>

bool try_adjacency(int a, int b, int pos_a_i, int pos_a_j, int pos_b_i, int pos_b_j, 
    int dist_border, int type_node_b, int GRID_SIZE, int TOTAL_GRID_SIZE, int NODE_SIZE,
    int t, pair<int,int> key, int global_pos_b, int *grid_place, int *pos_i, int *pos_j, 
    int *node_degree, int *grid_freedom, int *pos_io_i, int *pos_io_j, 
    int IN_OUT_SIZE, map<pair<int,int>,int> &edges_cost_local, int &cost_place, int arch, int rand) {

    if (rand)
        random_shuffle(ADJACENCY_FIRST.begin(), ADJACENCY_FIRST.end());
    
    int pos_b, cost;           
    for(int j = 0; j < ADJACENCY_SIZE; ++j) {
        if (j < 4) {
            pos_b_i = pos_a_i + ADJACENCY_FIRST[j].first;
            pos_b_j = pos_a_j + ADJACENCY_FIRST[j].second;
        } else {
            pos_b_i = pos_a_i + ADJACENCY[j][0];
            pos_b_j = pos_a_j + ADJACENCY[j][1];
        }

        if (pos_b_i < 0 || pos_b_j < 0 || pos_b_i >= GRID_SIZE || pos_b_j >= GRID_SIZE)
            continue;
        
        if (dist_border > 0) {
            cost = cost_local(pos_a_i, pos_a_j, pos_b_i, pos_b_j, GRID_SIZE, arch);
            if (cost > dist_border) continue;
        }

        pos_b = pos_b_i * GRID_SIZE + pos_b_j;

        if (placed(type_node_b, pos_b_i, pos_b_j, pos_b, GRID_SIZE, TOTAL_GRID_SIZE, 
            global_pos_b, grid_place, pos_i, pos_j, node_degree, grid_freedom)) {

            cost = cost_local(pos_a_i, pos_a_j, pos_b_i, pos_b_j, GRID_SIZE, arch);

            edges_cost_local[key] = cost;
            cost_place += cost;
            return true;
        }  
    }

    return false;
}

#endif
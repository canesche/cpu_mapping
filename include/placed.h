#ifndef __PLACED_H
#define __PLACED_H

#include <vector>

void create_freedrom_degree(const int GRID_SIZE, int *grid_freedom) {

    int i, j, count;

    #if __ARCH == 0
        for(i = 0; i < GRID_SIZE; ++i) {
            for(j = 0; j < GRID_SIZE; ++j) {
                count = 0;
                if (i-1 >= 0) ++count;
                if (i+1 < GRID_SIZE) ++count;
                if (j-1 >= 0) ++count;
                if (j+1 < GRID_SIZE) ++count;
                grid_freedom[i*GRID_SIZE+j] = count;
            }
        }
    #elif __ARCH == 1
        for(i = 0; i < GRID_SIZE; ++i) {
            for(j = 0; j < GRID_SIZE; ++j) {
                count = 0;
                if (i-1 >= 0) ++count;
                if (i-2 >= 0) ++count;
                if (i+1 < GRID_SIZE) ++count;
                if (i+2 < GRID_SIZE) ++count;
                if (j-1 >= 0) ++count;
                if (j-2 >= 0) ++count;
                if (j+1 < GRID_SIZE) ++count;
                if (j+2 < GRID_SIZE) ++count;
                grid_freedom[i*GRID_SIZE+j] = count;
            }
        }
    #else
        printf("Architecture not defined!\n");
    #endif
}

void update_grid_freedom(int pos, int *grid_freedom, const int GRID_SIZE, 
    const int TOTAL_GRID_SIZE) {
    
    int pos_i = pos / GRID_SIZE;
    int pos_j = pos % GRID_SIZE;

    int degree, x, y, aux_pos;
    grid_freedom[pos] = 0;

    for (int i = 0; i < SIZE_POS_TIPS_AUX; ++i){
        x = pos_i + POS_TIPS_AUX[i][0];
        y = pos_j + POS_TIPS_AUX[i][1];
        aux_pos = x*GRID_SIZE + y;
        if (x > 0 && x < GRID_SIZE && y > 0 && y < GRID_SIZE && grid_freedom[aux_pos] > 0) {
            degree = grid_freedom[aux_pos];
            if (degree > 1) grid_freedom[aux_pos] = degree-1;
        }
    }
}

void update_grid_freedom(int pos, int *grid_freedom, const int GRID_SIZE, 
    const int TOTAL_GRID_SIZE, vector<pair<int,int>> &aux) {
    
    aux.resize(0);
    int pos_i = pos / GRID_SIZE;
    int pos_j = pos % GRID_SIZE;

    int degree, x, y, aux_pos;
    grid_freedom[pos] = 0;

    for (int i = 0; i < SIZE_POS_TIPS_AUX; ++i){
        x = pos_i + POS_TIPS_AUX[i][0];
        y = pos_j + POS_TIPS_AUX[i][1];
        aux_pos = x*GRID_SIZE + y;
        if (x > 0 && x < GRID_SIZE && y > 0 && y < GRID_SIZE && grid_freedom[aux_pos] > 0) {
            degree = grid_freedom[aux_pos];
            if (degree > 1){ 
                grid_freedom[aux_pos] = degree-1;
                aux.push_back({degree-1,aux_pos});
            }
        }
    }
}

bool placed(int type_node, int pos_node_i, int pos_node_j, int pos_node, int GRID_SIZE, 
    int TOTAL_GRID_SIZE, int global_pos, int *grid_place, int *pos_i, int *pos_j, 
    int *node_degree, int *grid_freedom) {
    
    if (type_node == 0) {
        if (grid_place[pos_node] == type_node) {
            grid_place[pos_node] = -1;
            pos_i[global_pos] = pos_node_i;
            pos_j[global_pos] = pos_node_j;
#if __NEIGHBOURHOOD > 0
            node_degree[global_pos] -= 1;
            update_grid_freedom(pos_node, grid_freedom, GRID_SIZE, TOTAL_GRID_SIZE);
#endif
            return true;
        }
    } else { // type_node == 1 (IO)
        for (int j = 0; j < __LAYERS; j++) {
            int local_grid = pos_node+j*TOTAL_GRID_SIZE;
            if (grid_place[local_grid] == type_node) {
                grid_place[local_grid] = -1;
                pos_i[global_pos] = pos_node_i;
                pos_j[global_pos] = pos_node_j;
#if __NEIGHBOURHOOD > 0
                node_degree[global_pos] -= 1;
                update_grid_freedom(pos_node, grid_freedom, GRID_SIZE, TOTAL_GRID_SIZE);
#endif
                return true;
            }
        }
    }
    return false;
}

#endif
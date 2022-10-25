#ifndef __TYPE_NODE__H
#define __TYPE_NODE__H

void get_node_degree(Graph g, const int times, const int NODE_SIZE, int *node_degree){
    for (int i = 0; i < NODE_SIZE; ++i) {
        node_degree[i] = g.get_fanin()[i].size()+g.get_fanout()[i].size();
    }
    for (int j = 1; j < times; ++j) {
        for (int i = 0; i < NODE_SIZE; ++i) {
            node_degree[i+j*NODE_SIZE] = node_degree[i];
        }
    }
}

void get_type_node(vector<int> inputs, vector<int> outputs, int* type_node, 
    const int NODE_SIZE, int* type_matrix, const int GRID_SIZE,
    vector<int>& pos_input_output) {

    for (int i = 0; i < NODE_SIZE; ++i) {
        bool find = false;
        for (int j = 0; j < inputs.size(); ++j) {
            if (i == inputs[j]) {
                type_node[i] = 1;
                find = true;
                break;
            }
        }
        if (find) continue;
        for (int j = 0; j < outputs.size(); ++j) {
            if (i == outputs[j]) {
                type_node[i] = 1;
                find = true;
                break;
            }
        }
        if (find) continue;
        type_node[i] = 0;
    }

    for (int i = 0; i < GRID_SIZE*GRID_SIZE; ++i) {
        if (i == 0 || i == GRID_SIZE-1 || i == GRID_SIZE*(GRID_SIZE-1) || i == GRID_SIZE*GRID_SIZE-1) {
            for (int k = 0; k < __LAYERS; ++k) {
                type_matrix[i+k*GRID_SIZE*GRID_SIZE] = -1;
            }
        } else if (i < GRID_SIZE || i > GRID_SIZE*(GRID_SIZE-1) || i%GRID_SIZE==0 || (i+1)%GRID_SIZE==0) { 
            for (int k = 0; k < __LAYERS; ++k) {
                type_matrix[i+k*GRID_SIZE*GRID_SIZE] = 1;
            }
            pos_input_output.push_back(i);
        }
        else { 
            for (int k = 0; k < __LAYERS; ++k) {
                type_matrix[i+k*GRID_SIZE*GRID_SIZE] = 0;
            }
        }
    }
}

#endif
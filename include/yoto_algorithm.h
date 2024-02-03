#ifndef YOTO_ALGORITHM_H
#define YOTO_ALGORITHM_H

#include <Graph.h>
#include <get_critical_path.h>
#include <stack>
#include <algorithm>
#include <vector>

#include <list_adjacency.h>
#include <try_adjacency.h>
#include <get_cost.h>

void yoto_algorithm(const int NODE_SIZE, const int GRID_SIZE, const int EDGE_SIZE, 
        const int start, const int times, vector<tuple3> &vector_edges, int *node_degree,
        int *grid_freedom_original, int *pos_i, int* pos_j, int *type_node, int *type_matrix, 
        int *pos_io_i, int *pos_io_j, const int IN_OUT_SIZE, 
        vector<map<pair<int,int>,int>> &edges_cost, int *vector_cost, 
        int *list_borders, int arch);

void routing_yott(
    int *N,
    int *M,
    int EDGE_SIZE, 
    int NODE_SIZE, 
    vector<tuple3> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    int *pos_i, 
    int *pos_j, 
    int times,
    const int ARCH,
    int &bad_route,
    bool *successfulRoutings,
    int* vector_cost
);

bool verify_size_input(int GRID_SIZE, int nodes);

void yoto_main(
    string path_dot,
    int times,
    const int arch,
    const int parallel_mode = 0
) {
    string name = "";
    vector<string> aux;
    int max_threads = 1;

    aux = split(path_dot, "/");
    name = split(aux[aux.size()-1], ".")[0];

    if (parallel_mode == 1) {
        max_threads = omp_get_max_threads(); //get number of threads
        times += (max_threads - times % max_threads) % max_threads;
    }

    time_point<high_resolution_clock> start_clock, end_clock;
    double time_placement = 0, time_list, time_cp = 0;

    // Constructor Graph
    Graph g(path_dot);

    const unsigned int NODE_SIZE = g.num_nodes();
    const unsigned int EDGE_SIZE = g.num_edges();
    const unsigned int GRID_SIZE = ceil(sqrt(NODE_SIZE-g.get_inputs().size()-g.get_outputs().size())) + 2;
    const unsigned int TOTAL_GRID_SIZE = GRID_SIZE * GRID_SIZE;

    int *N = new int[times];
    int *M = new int[times];

    for (int i = 0; i < times; ++i) {
        N[i] = GRID_SIZE; //ceil(sqrt(nodes));
        M[i] = GRID_SIZE; //ceil(sqrt(nodes));
    }

    if (!verify_size_input(GRID_SIZE, g.get_inputs().size()+g.get_outputs().size())){
        printf("%s, Not solution for this GRID_SIZE %d\n", name.c_str(), GRID_SIZE);
        return;
    }

    vector<pair<int,int>> edges = g.get_edges();

    int* type_node = new int[NODE_SIZE];
    int* type_matrix = new int[GRID_SIZE*GRID_SIZE];

    vector<int> aux_in_out;

    get_type_node(g.get_inputs(), g.get_outputs(), type_node, NODE_SIZE, type_matrix, GRID_SIZE, aux_in_out);

    const int IN_OUT_SIZE = aux_in_out.size();
    int pos_io_i[IN_OUT_SIZE], pos_io_j[IN_OUT_SIZE];
    for (int i = 0; i < IN_OUT_SIZE; ++i) { 
        pos_io_i[i] = aux_in_out[i] / GRID_SIZE;
        pos_io_j[i] = aux_in_out[i] % GRID_SIZE;
    }

    vector<tuple3> vector_edges;
    vector<map<pair<int, int>,vector<tuple3>>> dic_CYCLE;
    /** generate start sequence:
        0: greedy neihborhood (large fanin/fanout)
        1: betweens centrality 
        2: Critical path
        3: random (ZigZag switch)
        4 or otherwise: zigzag without help preprocessing
    **/
    int change = 4;
    int *node_degree = new int[times*NODE_SIZE];
    int *list_borders = new int[NODE_SIZE];

    create_list_borders(g, NODE_SIZE, GRID_SIZE, list_borders);

    //for (int i = 0; i < NODE_SIZE; ++i) printf("%d %d\n", i, list_borders[i]);

    start_clock = high_resolution_clock::now();

    smart_transversal_algorithm(g, edges, NODE_SIZE, vector_edges, dic_CYCLE, change, times);

    end_clock = high_resolution_clock::now();
    time_list = duration_cast<Milliseconds>(end_clock - start_clock).count();

    int cost = 99999999;

    int *grid_freedom;
    int *pos_i = new int[times*NODE_SIZE];
    int *pos_j = new int[times*NODE_SIZE];

    for(int i = 0; i < times*NODE_SIZE; ++i) 
        pos_i[i] = pos_j[i] = -1;

    vector<map<pair<int,int>,int>> edges_cost;
    int *vector_cost = new int[times];

    if (parallel_mode == 0) {
        grid_freedom = new int[TOTAL_GRID_SIZE];

        create_freedrom_degree(GRID_SIZE, grid_freedom);

        start_clock = high_resolution_clock::now();
        yoto_algorithm(NODE_SIZE, GRID_SIZE, EDGE_SIZE, 0, times, vector_edges, 
            node_degree, grid_freedom, pos_i, pos_j, type_node, type_matrix, pos_io_i,
            pos_io_j, IN_OUT_SIZE, edges_cost, vector_cost, list_borders, arch);
        end_clock = high_resolution_clock::now();

        time_placement = duration_cast<Milliseconds>(end_clock-start_clock).count();
    } else {
        //Init arrays of each thread
        int **grid_freedom_thread = new int*[max_threads];
        vector<map<pair<int,int>,int>> *edges_cost_threads = new vector<map<pair<int,int>,int>>[max_threads];

        for(int i = 0; i < max_threads; i++) {
            grid_freedom_thread[i] = new int[TOTAL_GRID_SIZE];
            create_freedrom_degree(GRID_SIZE, grid_freedom_thread[i]);
        }

        int block_threads = times / max_threads;

        start_clock = high_resolution_clock::now();
        #pragma omp parallel for
        for(int i = 0; i < max_threads; i++){
            yoto_algorithm(NODE_SIZE, GRID_SIZE, EDGE_SIZE, i*block_threads,
            times/max_threads + (i < times % max_threads), vector_edges, 
            node_degree, grid_freedom_thread[i], pos_i, pos_j, type_node, type_matrix, 
            pos_io_i, pos_io_j, IN_OUT_SIZE, edges_cost_threads[i], vector_cost, 
            list_borders, arch);
        }
        end_clock = high_resolution_clock::now();
        time_placement = duration_cast<Milliseconds>(end_clock - start_clock).count();

        for(int i = 0; i < max_threads; i++){
            edges_cost.insert(edges_cost.end(), edges_cost_threads[i].begin(), edges_cost_threads[i].end()); 
        }
    }

    int** grid = new int*[times];
    for (int i = 0; i < times; ++i) {
        grid[i] = new int[TOTAL_GRID_SIZE];
    }

    for (int i = 0; i < times; ++i) {
        for (int j = 0; j < TOTAL_GRID_SIZE; ++j) {
            grid[i][j] = -1;
        }
        for (int j = 0; j < NODE_SIZE; ++j) {
            int pos_local = i*NODE_SIZE+j;
            int pos_node = pos_i[pos_local] * GRID_SIZE + pos_j[pos_local]; 
            grid[i][pos_node] = j;
        }
    }

    bool *successfullRoutings = new bool[times]; 
    for (int i = 0; i < times; ++i) {
        successfullRoutings[i] = true;
    }

    int bad_route = 0;

    // routing function
    routing_yott(N, M, EDGE_SIZE, NODE_SIZE, vector_edges, edges_cost, pos_i, pos_j, times, 0, bad_route, successfullRoutings, vector_cost);

    statistic(g, name, times, type_node, vector_edges, edges_cost, vector_cost);

    printf("%.2lf, %.2lf\n",time_placement,time_list);

    int best_index = get_best_index_yott(times, EDGE_SIZE, vector_edges, edges_cost, successfullRoutings);

    printGrid(grid[best_index], N[best_index], M[best_index]);

    output_graph_yott(vector_edges, edges_cost, best_index, EDGE_SIZE, g, name, "yoto", times, arch);

    // print the edges values and names
    for (auto e : edges) {
        printf("%d (%s) -> %d (%s)\n", e.first, g.get_name_node(e.first).c_str(), e.second, g.get_name_node(e.second).c_str());
    }
}

void yoto_algorithm(
    const int NODE_SIZE, 
    const int GRID_SIZE, 
    const int EDGE_SIZE, 
    const int start, 
    const int times, 
    vector<tuple3> &vector_edges, 
    int *node_degree,
    int *grid_freedom_original, 
    int *pos_i, 
    int* pos_j, 
    int *type_node, 
    int *type_matrix, 
    int *pos_io_i, 
    int *pos_io_j, 
    const int IN_OUT_SIZE, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    int *vector_cost, 
    int *list_borders,
    int arch
){
    
    int cost = -1, dist_border;
    const int TOTAL_GRID_SIZE = GRID_SIZE * GRID_SIZE;
    int pos_tips_size = 0, pos_inter_size = 0, pos_aux_size = 0, pos_forward_size = 0;
    pair<int, int> key, next_key;
    bool found;
    int a, b, pos_node, dist_max, degree, pos_a_i, pos_a_j, pos_b_i, pos_b_j;
    pair<int,int> biggest;
    vector<pair<int,int>> pos_tips, pos_forward, pos_inter, pos_aux;
    vector<int> pos_cost, pos_node_try;
    int size_pos_node_try, count, cost_min, cost_place, global_pos_a, global_pos_b, local_edge;
    int best_cost = 9999999, local_grid, rand_pos, pos_a, pos_b;
    int* grid_place = new int[TOTAL_GRID_SIZE];
    int* grid_freedom = new int[TOTAL_GRID_SIZE];
    int type_node_a, type_node_b, pos_cycle_size;
    map<pair<int,int>,int> edges_cost_local;
    vector<tuple3> pos_cycle;

    for (int t = 0; t < times; ++t) {
        pos_a = pos_b = -1;

        edges_cost_local.clear();

        for(int i = 0; i < TOTAL_GRID_SIZE; ++i) grid_place[i] = type_matrix[i];

        cost_place = 0;
        for (int i = 0; i < EDGE_SIZE; ++i) {
            
            local_edge = i+(t+start)*EDGE_SIZE;

            a = vector_edges[local_edge].v0;
            b = vector_edges[local_edge].v1;

            type_node_a = type_node[a];
            type_node_b = type_node[b];

            global_pos_a = a+(t+start)*NODE_SIZE;
            global_pos_b = b+(t+start)*NODE_SIZE;
            
            pos_a_i = pos_i[global_pos_a];
            pos_a_j = pos_j[global_pos_a];
            pos_b_i = pos_i[global_pos_b];
            pos_b_j = pos_j[global_pos_b];

            pos_a = pos_a_i * GRID_SIZE + pos_a_j;
            pos_b = pos_b_i * GRID_SIZE + pos_b_j;

            key = make_pair(a, b);
            
            if (pos_a_i != -1 && pos_b_i != -1) {
                cost = cost_local(pos_a_i, pos_a_j, pos_b_i, pos_b_j, GRID_SIZE, arch);
                cost_place += cost;
                edges_cost_local[key] = cost;
                continue;
            }
            
            size_pos_node_try = 0;

            // verify if side 'a' already place 
            if (pos_a_i == -1 && pos_b_i == -1) {
                // get the position random
                while (true) {
                    rand_pos = rand() % (IN_OUT_SIZE);
                    pos_a_i = pos_io_i[rand_pos];
                    pos_a_j = pos_io_j[rand_pos];
                    pos_a = pos_a_i * GRID_SIZE + pos_a_j;
                    //printf("pos_a (%d,%d)\n", pos_a_i, pos_a_j);
                    if (placed(type_node_a, pos_a_i, pos_a_j, pos_a, GRID_SIZE, TOTAL_GRID_SIZE, 
                        global_pos_a, grid_place, pos_i, pos_j, node_degree, grid_freedom)) {
                        break;
                    }
                }
            } 

            if (pos_b_i == -1) { // verify if side 'b' already place
                // try more closed
                if (try_adjacency(a, b, pos_a_i, pos_a_j, pos_b_i, pos_b_j, dist_border, 
                type_node_b, GRID_SIZE, TOTAL_GRID_SIZE, NODE_SIZE, t, key, global_pos_b, 
                grid_place, pos_i, pos_j, node_degree, grid_freedom, pos_io_i, pos_io_j, 
                IN_OUT_SIZE, edges_cost_local, cost_place, arch, 0)) continue;
            } 

            if (pos_a == -1 || pos_b == -1) {
                cost_place += 9999999;
                edges_cost_local[key] = 1000;
            }
        }
        edges_cost.push_back(edges_cost_local);
        vector_cost[t+start] = cost_place;
    }
    return;
}

#endif
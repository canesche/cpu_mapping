#ifndef YOTT_ALGORITHM_H
#define YOTT_ALGORITHM_H

#include <Graph.h>
#include <get_critical_path.h>
#include <stack>
#include <algorithm>
#include <vector>

#include <list_adjacency.h>
#include <try_adjacency.h>
#include <get_cost.h>

#define __LAYERS 1

void yott_algorithm(const int NODE_SIZE, const int GRID_SIZE, const int EDGE_SIZE, 
        const int start, const int times, vector<tuple3> &vector_edges, int *node_degree,
        int *grid_freedom_original, int *pos_i, int* pos_j, int *type_node, int *type_matrix, 
        int *pos_io_i, int *pos_io_j, const int IN_OUT_SIZE, 
        vector<map<pair<int,int>,int>> &edges_cost, int *vector_cost, 
        int *list_borders, vector<map<pair<int, int>,vector<tuple3>>> &dic_CYCLE);

void smart_transversal_algorithm(Graph g, vector<pair<int,int>> edges, 
    const int NODE_SIZE, vector<tuple3> &EDGES, vector<map<pair<int, int>,vector<tuple3>>> &dic_CYCLE, 
    int change, const int times);

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
    bool *successfulRoutings
);

bool verify_size_input(int GRID_SIZE, int nodes);

void yott_main(
    string path_dot,
    int times,
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
    int* type_matrix = new int[GRID_SIZE*GRID_SIZE*__LAYERS];

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
    const int change = 5;
    int *node_degree = new int[times*NODE_SIZE];
    int *list_borders = new int[NODE_SIZE];

    //g.print_graph_number();
    #if __NEIGHBOURHOOD > 0
    get_node_degree(g, times, NODE_SIZE, node_degree);
    #endif

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
        yott_algorithm(NODE_SIZE, GRID_SIZE, EDGE_SIZE, 0, times, vector_edges, 
            node_degree, grid_freedom, pos_i, pos_j, type_node, type_matrix, pos_io_i,
            pos_io_j, IN_OUT_SIZE, edges_cost, vector_cost, list_borders, dic_CYCLE);
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
            yott_algorithm(NODE_SIZE, GRID_SIZE, EDGE_SIZE, i*block_threads,
            times/max_threads + (i < times % max_threads), vector_edges, 
            node_degree, grid_freedom_thread[i], pos_i, pos_j, type_node, type_matrix, 
            pos_io_i, pos_io_j, IN_OUT_SIZE, edges_cost_threads[i], vector_cost, 
            list_borders, dic_CYCLE);
        }
        end_clock = high_resolution_clock::now();
        time_placement = duration_cast<Milliseconds>(end_clock - start_clock).count();

        for(int i = 0; i < max_threads; i++){
            edges_cost.insert(edges_cost.end(), edges_cost_threads[i].begin(), edges_cost_threads[i].end()); 
        }
    }

    bool *successfullRoutings = new bool[times]; 
    for (int i = 0; i < times; ++i) {
        successfullRoutings[i] = true;
    }

    
    int bad_route = 0;

    // routing function
    routing_yott(N, M, EDGE_SIZE, NODE_SIZE, vector_edges, edges_cost, pos_i, pos_j, times, 0, bad_route, successfullRoutings);

    statistic(g, name, times, type_node, vector_edges, edges_cost, vector_cost);

    //printf("\n%d\n", best_idx);

    printf("%.2lf,%.2lf\n",time_placement,time_list);

    int best_index = get_best_index_yott(times, EDGE_SIZE, vector_edges, edges_cost, successfullRoutings);

    output_graph_yott(vector_edges, edges_cost, best_index, EDGE_SIZE, g);
}

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
    bool *successfulRoutings
) {

    int a, b;
    int pos_a_x, pos_a_y, pos_b_x, pos_b_y;
    int diff_x, diff_y, dist_x, dist_y;
    int pos_node_x, pos_node_y;
    int change, pe_curr, step;
    int GRID_X, GRID_Y, TOTAL_GRID_SIZE;
    bool routing;
    int dFreedom = 0;

    bool neighbor_right, neighbor_left, neighbor_down, neighbor_top;
    bool neighbor_right2, neighbor_left2, neighbor_down2, neighbor_top2;

    if (ARCH == 0) // MESH
        dFreedom = 4;
    else if (ARCH == 1 || ARCH == 2) // 1-HOP OR CHESS
        dFreedom = 8;
    else {
        printf("Architecture not defined!\n");
        exit(1);    
    }
    
    int count_per_curr[TOTAL_GRID_SIZE];
    int grid[TOTAL_GRID_SIZE][dFreedom];

    for (int t = 0; t < times; ++t) {

        routing = true;
        successfulRoutings[t] = true;
        GRID_X = N[t];
        GRID_Y = M[t];

        TOTAL_GRID_SIZE = GRID_X * GRID_Y;
        vector<vector<map<int,int>>> grid_route(TOTAL_GRID_SIZE, vector<map<int,int>>(dFreedom));
        int count_per_curr[TOTAL_GRID_SIZE];

        vector<edge_route> edge(EDGE_SIZE);
        vector<vector<int>> path(NODE_SIZE);

        for (int i = 0; i < EDGE_SIZE; ++i) {

            a = vector_edges[t*EDGE_SIZE+i].v0;
            b = vector_edges[t*EDGE_SIZE+i].v1;
            pos_a_x = pos_i[t*NODE_SIZE+a];
            pos_a_y = pos_j[t*NODE_SIZE+a];
            pos_b_x = pos_i[t*NODE_SIZE+b];
            pos_b_y = pos_j[t*NODE_SIZE+b];

            //printf("%2d (%2d,%2d) -> %2d (%2d,%2d)\n", a, pos_a_x, pos_a_y, b, pos_b_x, pos_b_y);

            diff_x = pos_b_x - pos_a_x;
            diff_y = pos_b_y - pos_a_y;

            dist_x = abs(diff_x);
            dist_y = abs(diff_y);

            pos_node_x = pos_a_x;
            pos_node_y = pos_a_y;

            change = 0;
            step = 0;

            for (int j = 0; j < TOTAL_GRID_SIZE; ++j)
                count_per_curr[j] = 0; 
            
            while(dist_x != 0 || dist_y != 0) {
                // print_grid_route(grid_route, GRID_X, GRID_Y, dFreedom);

                change = 0;
                diff_x = pos_b_x - pos_node_x;
                diff_y = pos_b_y - pos_node_y;

                ++step; // one step

                //printf("%d %d\n", diff_x, diff_y);
                  
                // get the current position node
                pe_curr = pos_node_x * GRID_Y + pos_node_y;
                count_per_curr[pe_curr]++;
                path[a].push_back(pe_curr);
                
                // [pe], [0 = top, 1 = right, 2 = down, 3 = left, 4 = top2, 5 = right2, 6 = down2, 7 = left2]
                     
                if (ARCH == 1 || ARCH == 2) {
                    neighbor_right2 = pe_curr+2 < (pos_node_x+1)*GRID_Y && (grid_route[pe_curr][5].count(step) == 0 || ((grid_route[pe_curr][5].count(step) > 0) ? grid_route[pe_curr][5][step] == a : true)) && ((grid_route[pe_curr+2][7].count(step) > 0)        ? grid_route[pe_curr+2][7][step] != a : true)        && count_per_curr[pe_curr+2] == 0;
                    neighbor_left2 = pe_curr-2 >= (pos_node_x)*GRID_Y   && (grid_route[pe_curr][7].count(step) == 0 || ((grid_route[pe_curr][7].count(step) > 0) ? grid_route[pe_curr][7][step] == a : true)) && ((grid_route[pe_curr-2][5].count(step) > 0)        ? grid_route[pe_curr-2][5][step] != a : true)        && count_per_curr[pe_curr-2] == 0;
                    neighbor_down2 = pe_curr+2*GRID_Y < TOTAL_GRID_SIZE && (grid_route[pe_curr][6].count(step) == 0 || ((grid_route[pe_curr][6].count(step) > 0) ? grid_route[pe_curr][6][step] == a : true)) && ((grid_route[pe_curr+2*GRID_Y][4].count(step) > 0) ? grid_route[pe_curr+2*GRID_Y][4][step] != a : true) && count_per_curr[pe_curr+2*GRID_Y] == 0;
                    neighbor_top2 = pe_curr-2*GRID_Y >= 0               && (grid_route[pe_curr][4].count(step) == 0 || ((grid_route[pe_curr][4].count(step) > 0) ? grid_route[pe_curr][4][step] == a : true)) && ((grid_route[pe_curr-2*GRID_Y][6].count(step) > 0) ? grid_route[pe_curr-2*GRID_Y][6][step] != a : true) && count_per_curr[pe_curr-2*GRID_Y] == 0;
                }

                neighbor_right = pe_curr+1 < (pos_node_x+1)*GRID_Y && (grid_route[pe_curr][1].count(step) == 0 || ((grid_route[pe_curr][1].count(step) > 0) ? grid_route[pe_curr][1][step] == a : true)) && ((grid_route[pe_curr+1][3].count(step) > 0)        ? grid_route[pe_curr+1][3][step] != a : true)        && count_per_curr[pe_curr+1] == 0;
                neighbor_left = pe_curr-1 >= (pos_node_x)*GRID_Y   && (grid_route[pe_curr][3].count(step) == 0 || ((grid_route[pe_curr][3].count(step) > 0) ? grid_route[pe_curr][3][step] == a : true)) && ((grid_route[pe_curr-1][1].count(step) > 0)        ? grid_route[pe_curr-1][1][step] != a : true)        && count_per_curr[pe_curr-1] == 0;
                neighbor_down = pe_curr+GRID_Y < TOTAL_GRID_SIZE   && (grid_route[pe_curr][2].count(step) == 0 || ((grid_route[pe_curr][2].count(step) > 0) ? grid_route[pe_curr][2][step] == a : true)) && ((grid_route[pe_curr+GRID_Y][0].count(step) > 0)   ? grid_route[pe_curr+GRID_Y][0][step] != a : true)   && count_per_curr[pe_curr+GRID_Y] == 0;
                neighbor_top = pe_curr-GRID_Y >= 0                 && (grid_route[pe_curr][0].count(step) == 0 || ((grid_route[pe_curr][0].count(step) > 0) ? grid_route[pe_curr][0][step] == a : true)) && ((grid_route[pe_curr-GRID_Y][2].count(step) > 0)   ? grid_route[pe_curr-GRID_Y][2][step] != a : true)   && count_per_curr[pe_curr-GRID_Y] == 0;

                //printf("r = %d l = %d d = %d t = %d\n", neighbor_right, neighbor_left, neighbor_down, neighbor_top);
                
                // go to right neighbor
                if (ARCH == 0 || (ARCH == 2 && pe_curr % 2 == 0)) { // ARCH MESH OR CHESS EVEN
                    if (diff_y > 0 && neighbor_right) {
                        grid_route[pe_curr][1][step] = a;
                        pos_node_y += 1;
                        change = 1;
                    } else if (diff_y < 0 && neighbor_left) {
                        grid_route[pe_curr][3][step] = a;
                        pos_node_y -= 1;
                        change = 1;
                    } else if (diff_x > 0 && neighbor_down) {
                        grid_route[pe_curr][2][step] = a;
                        pos_node_x += 1;
                        change = 1;
                    } else if (diff_x < 0 && neighbor_top) {
                        grid_route[pe_curr][0][step] = a;
                        pos_node_x -= 1;
                        change = 1;
                    }
                    
                    if (change == 0) { //tenta fazer uma voltinha e aumentar o caminho
                        if (neighbor_right) {
                            grid_route[pe_curr][1][step] = a;
                            pos_node_y += 1;
                            change = 1;
                        } 
                        else if (neighbor_left) {
                            grid_route[pe_curr][3][step] = a;
                            pos_node_y -= 1;
                            change = 1;
                        }
                        else if (neighbor_down) {
                            grid_route[pe_curr][2][step] = a;
                            pos_node_x += 1;
                            change = 1;
                        } 
                        else if (neighbor_top) {
                            grid_route[pe_curr][0][step] = a;
                            pos_node_x -= 1;
                            change = 1;
                        }
                    }
                } else if (ARCH == 1 || (ARCH == 2 && pe_curr % 2 != 0)) { // ARCH 1-HOP OR CHESS ODD
                    if (diff_y > 1 && neighbor_right2) {
                        grid_route[pe_curr][5][step] = a;
                        pos_node_y += 2;
                        change = 1;
                    } else if (diff_y > 0 && neighbor_right) {
                        grid_route[pe_curr][1][step] = a;
                        pos_node_y += 1;
                        change = 1;
                    } else if (diff_y < -1 && neighbor_left2) {
                        grid_route[pe_curr][7][step] = a;
                        pos_node_y -= 2;
                        change = 1;
                    } else if (diff_y < 0 && neighbor_left) {
                        grid_route[pe_curr][3][step] = a;
                        pos_node_y -= 1;
                        change = 1;
                    } else if (diff_x > 1 && neighbor_down2) {
                        grid_route[pe_curr][6][step] = a;
                        pos_node_x += 2;
                        change = 1;
                    } else if (diff_x > 0 && neighbor_down) {
                        grid_route[pe_curr][2][step] = a;
                        pos_node_x += 1;
                        change = 1;
                    } else if (diff_x < -1 && neighbor_top2) {
                        grid_route[pe_curr][4][step] = a;
                        pos_node_x -= 2;
                        change = 1;
                    } else if (diff_x < 0 && neighbor_top) {
                        grid_route[pe_curr][0][step] = a;
                        pos_node_x -= 1;
                        change = 1;
                    }

                    if (change == 0) { //tenta fazer uma voltinha e aumentar o caminho
                        if (neighbor_right) {
                            grid_route[pe_curr][1][step] = a;
                            pos_node_y += 1;
                            change = 1;
                        } 
                        else if (neighbor_left) {
                            grid_route[pe_curr][3][step] = a;
                            pos_node_y -= 1;
                            change = 1;
                        }
                        else if (neighbor_down) {
                            grid_route[pe_curr][2][step] = a;
                            pos_node_x += 1;
                            change = 1;
                        } 
                        else if (neighbor_top) {
                            grid_route[pe_curr][0][step] = a;
                            pos_node_x -= 1;
                            change = 1;
                        } else if (neighbor_right2) {
                            grid_route[pe_curr][5][step] = a;
                            pos_node_y += 2;
                            change = 1;
                        }
                        else if (neighbor_left2) {
                            grid_route[pe_curr][7][step] = a;
                            pos_node_y -= 2;
                            change = 1;
                        } 
                        else if (neighbor_down2) {
                            grid_route[pe_curr][6][step] = a;
                            pos_node_x += 2;
                            change = 1;
                        }
                        else if (neighbor_top2) {
                            grid_route[pe_curr][4][step] = a;
                            pos_node_x -= 2;
                            change = 1;
                        }
                    }
                }

                if (change == 0) { // not routing possible
                    bad_route++; 
                    routing = false;
                    step = 1000;
                    break;
                }

                dist_x = abs(pos_b_x - pos_node_x); 
                dist_y = abs(pos_b_y - pos_node_y);
            }

            if (!routing) { 
                successfulRoutings[t] = false;
                break;
            }

            // printf("%d -> %d old Cost: %d new Cost: %d\n", a, b, edge[i].cost, step);
            edges_cost[t][make_pair(a,b)] = step;
        }   
    }
}

bool verify_size_input(int GRID_SIZE, int nodes) {
    int sum = 0;
    for(int i = 0; i < GRID_SIZE*GRID_SIZE; ++i) {
        if (i == 0 || i == GRID_SIZE-1 || i == GRID_SIZE*(GRID_SIZE-1) || i == GRID_SIZE*GRID_SIZE-1) continue;
        if (i < GRID_SIZE || i > GRID_SIZE*(GRID_SIZE-1) || i%GRID_SIZE==0 || (i+1)%GRID_SIZE==0) sum += __LAYERS;
    }
    return sum >= nodes;
}

void get_tips(int pos_i, int pos_j, int *grid_freedom, const int GRID_SIZE, int node_type, 
    int *grid_place, int dist_border, vector<pair<int,int>> &aux, int *pos_io_i, 
    int *pos_io_j, int IN_OUT_SIZE) {
    
    aux.resize(0);
    
    int aux_pos, x, y, cost;

    if (dist_border == 0) { // node not closed to IO
        for (int i = 0; i < SIZE_POS_TIPS_AUX; ++i){
            x = pos_i + POS_TIPS_AUX[i][0];
            y = pos_j + POS_TIPS_AUX[i][1];
            aux_pos = x*GRID_SIZE + y;
            if (x > 0 && x < GRID_SIZE && y > 0 && y < GRID_SIZE 
#if __NEIGHBOURHOOD > 0
                && grid_freedom[aux_pos] > 0 
#endif                
                && grid_place[aux_pos] == node_type) {
                
                aux.push_back({grid_freedom[aux_pos],aux_pos});
            }
        }
    } else { // node close to IO
#if __DEBUG == 1
    printf("Pos Border: %d\n", dist_border);
#endif
        for (int i = 0; i < SIZE_POS_TIPS_AUX; ++i){
            x = pos_i + POS_TIPS_AUX[i][0];
            y = pos_j + POS_TIPS_AUX[i][1];
            aux_pos = x*GRID_SIZE + y;
            if (x > 0 && x < GRID_SIZE && y > 0 && y < GRID_SIZE 
#if __NEIGHBOURHOOD > 0
                && grid_freedom[aux_pos] > 0 
#endif
                && grid_place[aux_pos] == node_type) {
                for (int j = 0; j < IN_OUT_SIZE; ++j) {
                    cost = cost_local(x, y, pos_io_i[j], pos_io_j[j], GRID_SIZE);
                    if (cost <= dist_border) {
                        aux.push_back({grid_freedom[aux_pos],aux_pos});
                        break;
                    }
                }
            }
        }
    }
}

void get_pos_neighbor(int pos_node, int dist_max, int *grid_freedom, 
    const int GRID_SIZE, const int TOTAL_GRID_SIZE, vector<pair<int,int>> &aux) {
    
    aux.resize(0);

    int pos_i, pos_j, new_pos_i, new_pos_j, pos_global, degree;
    pos_i = pos_node / GRID_SIZE;
    pos_j = pos_node % GRID_SIZE;
    int diff;
    bool can;

    #if __ARCH == 1
        dist_max *= 2;
    #endif
    if (dist_max == 0) dist_max = 1;

    for (int i = -dist_max; i <= dist_max; ++i){
        for (int j = -dist_max; j <= dist_max; ++j) {
            can = false;
            new_pos_i = i + pos_i;
            new_pos_j = j + pos_j;
            
            #if __ARCH == 0
                diff = abs(i) + abs(j) - 1;
            #elif __ARCH == 1
                diff = (abs(i)/2 + abs(i)%2 + abs(j)/2 + abs(j)%2 - 1);
            #endif
            
            can = (diff >= 0 && diff < dist_max);

            pos_global = new_pos_i*GRID_SIZE + new_pos_j;
            degree = grid_freedom[pos_global];

            if (can && new_pos_i >= 0 && new_pos_j >= 0 &&
                new_pos_i < GRID_SIZE && new_pos_j < GRID_SIZE && degree > 0){
                aux.push_back({degree, pos_global});
            }
        }
    }
}

void intersection(vector<pair<int,int>> &pos_inter, vector<pair<int,int>> &pos_aux) {
    vector<pair<int,int>> pos_aux_inter;

    sort(pos_inter.begin(), pos_inter.end());
    sort(pos_aux.begin(), pos_aux.end());

    set_intersection(pos_inter.begin(), pos_inter.end(), pos_aux.begin(), pos_aux.end(), back_inserter(pos_aux_inter));

    // update pos_inter
    pos_inter = pos_aux_inter;
}

void subtraction(vector<pair<int,int>> &pos, vector<pair<int,int>> &pos_aux) {
    vector<pair<int,int>> pos_aux_sub;

    sort(pos.begin(), pos.end());
    sort(pos_aux.begin(), pos_aux.end());

    set_difference(pos.begin(), pos.end(), pos_aux.begin(), pos_aux.end(), back_inserter(pos_aux_sub));

    // update position
    pos = pos_aux_sub;
}

int best_place_degree(vector<pair<int,int>> &pos_aux, int pos_size, 
    int *grid_place, int pos_degree, int type_node_pos, int *type_node) {

    // pos_aux = pair<degree, pos>
    // degree fanin e fanout do no em questao
    int best_degree = pos_aux[0].first;
    int best_pos = pos_aux[0].second;
    int pos, degree;

#if __NEIGHBOURHOOD == 1 // greedy
    for (int j = 0; j < pos_size; ++j) {
        pos = pos_aux[j].second;
        degree = pos_aux[j].first;
        // condicao 1
        if (type_node_pos == grid_place[pos] && degree >= pos_degree) return pos;
    }
#elif __NEIGHBOURHOOD == 2 //search 
    for (int j = 0; j < pos_size; ++j) {
        pos = pos_aux[j].second;
        degree = pos_aux[j].first;
        // condicao 2
        //pos_degree = |PE|, degree = |b|
        if(type_node_pos == grid_place[pos] && best_degree >= degree && pos_degree >= degree) {//best_degree >= degree porque senão nunca considera o primeiro
            best_degree = degree;
            best_pos = pos;
        } 
    }
    return best_pos;
#elif __NEIGHBOURHOOD == 3 //search - limiar 1
    for (int j = 1; j < pos_size; ++j) {
        pos = pos_aux[j].second;
        degree = pos_aux[j].first;
        // condicao 3
        if(type_node_pos == grid_place[pos] && best_degree >= degree && pos_degree > degree) {
            best_degree = degree;
            best_pos = pos;
        } 
    }
    return best_pos;
#elif __NEIGHBOURHOOD >= 4
     for (int j = 1; j < pos_size; ++j) {
        pos = pos_aux[j].second;
        degree = pos_aux[j].first;

        if (pos_degree > 3) { // get the most degree 
            if (pos > best_degree && type_node_pos == grid_place[pos]) {
                best_degree = degree;
                best_pos = pos;
            }
        } else { 
            // get the minimun degree acceptable
            if (best_degree > degree && pos_degree >= degree && type_node_pos == grid_place[pos]) {
                best_degree = degree;
                best_pos = pos;
            }
        }
    }
#endif
    return best_pos;
}

void yott_algorithm(
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
    vector<map<pair<int, int>,vector<tuple3>>> &dic_CYCLE
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
    int* grid_place = new int[TOTAL_GRID_SIZE*__LAYERS];
    int* grid_freedom = new int[TOTAL_GRID_SIZE];
    int type_node_a, type_node_b, pos_cycle_size;
    map<pair<int,int>,int> edges_cost_local;
    vector<tuple3> pos_cycle;

    for (int t = 0; t < times; ++t) {
        pos_a = pos_b = -1;

#if __DEBUG == 1
        printf("---> try: %d\n",t);
#endif

        edges_cost_local.clear();
#if __NEIGHBOURHOOD > 0
        for(int i = 0; i < TOTAL_GRID_SIZE; ++i) grid_freedom[i] = grid_freedom_original[i];
#endif
        for(int i = 0; i < TOTAL_GRID_SIZE*__LAYERS; ++i) grid_place[i] = type_matrix[i];

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

#if __DEBUG == 1
            printf("\n---Operation: %d(%d,%d) -> %d (%d,%d)\n", a, pos_a_i, pos_a_j, b, pos_b_i, pos_b_j);
#endif
            
            if (pos_a_i != -1 && pos_b_i != -1) {
                cost = cost_local(pos_a_i, pos_a_j, pos_b_i, pos_b_j, GRID_SIZE);
                cost_place += cost;
                edges_cost_local[key] = cost;
#if __DEBUG == 1
                printf("PLACE TYPE ALREADY: %2d -> %2d Cost: %d\n", a, b, cost);
                print_grid_elements(pos_i, pos_j, t, NODE_SIZE, GRID_SIZE, a, b);
#endif
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

                        get_tips(pos_a_i, pos_a_j, grid_freedom, GRID_SIZE, type_node_b, grid_place, 0, pos_tips, pos_io_i, pos_io_j, IN_OUT_SIZE);
                        break;
                    }
                }
            } else { // get tips of pos_a
                get_tips(pos_a_i, pos_a_j, grid_freedom, GRID_SIZE, type_node_b, grid_place, list_borders[b], pos_tips, pos_io_i, pos_io_j, IN_OUT_SIZE);
            }

            if (pos_b_i == -1) { // verify if side 'b' already place
                found = false;
#if __NEIGHBOURHOOD > 0
                pos_tips_size = pos_tips.size();
#else
                pos_tips_size = 0;
#endif
#if __CYCLE == 1
                pos_cycle = dic_CYCLE[t][key];
                pos_cycle_size = pos_cycle.size();
#else
                pos_cycle_size = 0;
#endif
                pos_forward = pos_tips;
                pos_forward_size = 0;
                pos_node_try.clear();
#if __DEBUG == 1
                cout << "POS TIPS: ";
                for (int j = 0; j < pos_tips_size; ++j) {
                    cout << pos_tips[j].second << " ";
                }
                cout << endl;
#endif

                if (pos_cycle_size > 0 && pos_tips_size > 0){
                    // forward
                    if (i < EDGE_SIZE-1) {
                        int next_a = vector_edges[i+1].v0;
                        int next_b = vector_edges[i+1].v1;
                        next_key = func_key(next_a, next_b);

                        for (int j = 0, n = dic_CYCLE[t][next_key].size(); j < n; ++j) {
                            if (dic_CYCLE[t][next_key][j].v1 == 0) {
#if __DEBUG == 1
                                cout << "FORWARD ";
                                cout << dic_CYCLE[t][next_key][j].v0 << "d" << dic_CYCLE[t][next_key][j].v1 << endl;
#endif
                                pos_node = pos_a; //pos[dic_CYCLE[next_key][j].first];
                                
                                /*cout << "NODE " << dic_CYCLE[next_key][j].first << " pos= " << pos_node << " pos: ";
                                //pos_aux = get_tips(pos_node, grid_freedom, GRID_SIZE, TOTAL_GRID_SIZE, 0);
                                
                                get_pos_neighbor(pos_node, 0, grid_freedom, GRID_SIZE, TOTAL_GRID_SIZE, pos_aux);
                                
                                for (int k = 0; k < pos_aux.size(); ++k) {
                                    cout << pos_aux[k].second << " ";
                                }
                                cout << endl;
                                pos_aux_size = pos_aux.size();
                                subtraction(pos_forward, pos_aux);
                                pos_forward_size = pos_forward.size();
                                
                                cout << "pos forward: ";
                                for (int k = 0; k < pos_forward.size(); ++k) {
                                    cout << pos_forward[k].second << " ";
                                }
                                cout << endl;
                                */
                            }
                        }
                        /*
                        if (pos_forward_size > 0) {
                            found = place_local(pos_forward, pos_forward_size, grid_place, b, node_degree, 
                                        pos, grid_freedom, GRID_SIZE, TOTAL_GRID_SIZE, arch, a, 5);
                            
                            cout << "place b = " << b << " pos= " << pos[b] << endl;
                            // skip to next edge
                            if (found) continue;
                        }*/
                    }

                    pos_inter = pos_forward;
                    pos_inter_size = pos_forward.size();

#if __DEBUG == 1                  
                    if (pos_cycle_size > 0) cout << "\nLIST CYCLE: ";
                    for (int j = 0; j < pos_cycle_size; ++j) cout << pos_cycle[j].v0 << "d" << pos_cycle[j].v1 << " ";
                    cout << endl;
#endif

                    for (int j = 0; j < pos_cycle_size; ++j) {
                        
                        if (pos_cycle[j].v1 > __MAX_ANNOTATION) continue;

                        pos_node_try.push_back(pos_cycle[j].v0);
                        
                        pos_node = pos_i[pos_cycle[j].v0] * GRID_SIZE + pos_j[pos_cycle[j].v0];
                        dist_max = pos_cycle[j].v1;
                        
                        get_pos_neighbor(pos_node, dist_max, grid_freedom, GRID_SIZE, TOTAL_GRID_SIZE, pos_aux);

#if DEBUG == 1
                        cout << "Node " << pos_cycle[j].v0 << "d" << dist_max << " pos= " << pos[pos_cycle[j].v0] << endl;
                        for (int k = 0; k < pos_aux.size(); ++k) {
                            cout << pos_aux[k].second << " ";
                        }
                        cout << endl;
#endif

                        pos_aux_size = pos_aux.size();
                        intersection(pos_inter, pos_aux);
                        pos_inter_size = pos_inter.size();

#if __DEBUG == 1
                        cout << "intersection: ";
                        for (int j = 0; j < pos_inter_size; ++j) {
                            cout << pos_inter[j].second << " "; 
                        }
                        cout << endl;
#endif
                        if (pos_inter_size == 0) {
#if __DEBUG == 1
                            cout << "VAZIO INTERSECAO" << endl;
#endif
                            break;
                        } 
                    }
                }

                if (pos_inter_size > 0) {
                    
                    pos_b = best_place_degree(pos_inter, pos_inter_size, grid_place, 
                        node_degree[global_pos_b], type_node_b, type_node);
                    pos_b_i = pos_b / GRID_SIZE;
                    pos_b_j = pos_b % GRID_SIZE;

                    if (placed(type_node_b, pos_b_i, pos_b_j, pos_b, GRID_SIZE, TOTAL_GRID_SIZE, 
                        global_pos_b, grid_place, pos_i, pos_j, node_degree, grid_freedom)) {
                        
                        cost = cost_local(pos_a_i, pos_a_j, pos_b_i, pos_b_j, GRID_SIZE);

                        #if __DEBUG == 1
                            printf("PLACE TYPE INTERCESSAO: %2d -> %2d Cost: %d\n", a, b, cost);
                            print_grid_elements(pos_i, pos_j, t, NODE_SIZE, GRID_SIZE, a, b);
                        #endif

                        cost_place += cost;
                        edges_cost_local[key] = cost;
                        continue; // skip to next edge
                    }
                }

                if (pos_tips_size > 0) {
                    // Trying the get cost better
                    pos_aux = pos_tips;
                    size_pos_node_try = pos_node_try.size();
                    pos_cost.clear();
                    
                    if (size_pos_node_try > 0) {
                        for (int k = 0; k < size_pos_node_try; ++k) {
                            for (int j = 0; j < pos_tips_size; ++j) {
                                int aux = pos_node_try[k]+(t+start)*NODE_SIZE;
                                cost = cost_local(pos_i[aux]/GRID_SIZE,pos_j[aux]%GRID_SIZE,pos_tips[j].second/GRID_SIZE,pos_tips[j].second%GRID_SIZE, GRID_SIZE);
                                if (k == 0) {
                                    pos_cost.push_back(cost);
                                } else {
                                    pos_cost[j] += cost;
                                }
                            }
                        }
                        // choose the min cost
                        cost_min = pos_cost[0];
                        pos_aux.clear();
                        for (int j = 1; j < pos_tips_size; ++j) {
                            cost = pos_cost[j];
#if __DEBUG == 1
                            cout << "Pos = " << pos_tips[j].second << " weight " << cost << endl;
#endif
                            if (cost_min > cost) {
                                cost_min = cost;
                                pos_aux.clear();
                                pos_aux.push_back({pos_tips[j].first, pos_tips[j].second});
                            } else if (cost == cost_min) {
                                pos_aux.push_back({pos_tips[j].first, pos_tips[j].second});
                            }
                        }
                    }

                    pos_aux_size = pos_aux.size();
#if __DEBUG == 1
                    printf("Place TIPS: ");
                    for (int j = 0; j < pos_aux_size; ++j) printf("%d ", pos_aux[j].second);
                    printf("\n");
#endif

                    int pos_b = best_place_degree(pos_aux, pos_aux_size, grid_place, 
                        node_degree[global_pos_b], type_node_b, type_node);

                    pos_b_i = pos_b / GRID_SIZE;
                    pos_b_j = pos_b % GRID_SIZE;

                    if (placed(type_node_b, pos_b_i, pos_b_j, pos_b, GRID_SIZE, TOTAL_GRID_SIZE, 
                        global_pos_b, grid_place, pos_i, pos_j, node_degree, grid_freedom)) {

                        cost = cost_local(pos_a_i, pos_a_j, pos_b_i, pos_b_j, GRID_SIZE);

                        #if __DEBUG == 1
                            printf("PLACE TYPE DICAS: %2d -> %2d Cost: %d\n", a, b, cost);
                            print_grid_elements(pos_i, pos_j, t, NODE_SIZE, GRID_SIZE, a, b);
                        #endif

                        cost_place += cost;
                        edges_cost_local[key] = cost;
                        continue; // skip to next edge
                    }                                     
                }

                // try more closed
                if (try_adjacency(a, b, pos_a_i, pos_a_j, pos_b_i, pos_b_j, dist_border, 
                type_node_b, GRID_SIZE, TOTAL_GRID_SIZE, NODE_SIZE, t, key, global_pos_b, 
                grid_place, pos_i, pos_j, node_degree, grid_freedom, pos_io_i, pos_io_j, 
                IN_OUT_SIZE, edges_cost_local, cost_place)) continue;
            } 

            if (pos_a == -1 || pos_b == -1){
                //cerr << "No solution found!" << endl; 
                cost_place += 9999999;
                edges_cost_local[key] = 1000;
            }
        }

        // get the best
        edges_cost.push_back(edges_cost_local);
        vector_cost[t+start] = cost_place;
    }
    return;
}

#endif
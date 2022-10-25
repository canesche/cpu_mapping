#ifndef __ROUTING_H
#define __ROUTING_H

using namespace std;

#include <vector>
#include <map>
#include <algorithm>

struct edge_route {
    int a, b, cost;
    edge_route(int i = 0, int j = 0, int c = 0) : a(i), b(j), cost(c) {}
};

void print_grid_route(vector<vector<map<int,int>>> &grid, int N, int M, int dFreedom) {
    printf("\n");
    for (int i = 0; i < N*M; ++i) {
        printf("PE %d [", i);
        for (int j = 0; j < dFreedom; ++j) {
            printf("%d: ", j);
            for (auto m : grid[i][j]) {
                printf("{%d:%d}", m.first, m.second);
            }
            if (j != dFreedom-1) printf(", ");
        }
        printf("]\n");
    }
    printf("\n");
}

bool compare_edge(edge_route &a, edge_route &b) {
    return a.cost > b.cost;
}

void routing(int* N, int* M, int EDGE_SIZE, int NODE_SIZE, vector<pair<int,int>> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, int times, bool* successfulRoutings, int **grid, int **positions,
    int &bad_route, int ARCH) {
    
    int *pos_i = new int[NODE_SIZE];
    int *pos_j = new int[NODE_SIZE];

    int a, b;
    int pos_a_x, pos_a_y, pos_b_x, pos_b_y;
    int diff_x, diff_y, dist_x, dist_y;
    int pos_node_x, pos_node_y;
    int change, pe_curr, step;
    int GRID_X, GRID_Y, TOTAL_GRID_SIZE;
    int *gridFlat;
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
    
    for (int t = 0; t < times; ++t) {
        
        routing = true;
        successfulRoutings[t] = true;
        GRID_X = N[t];
        GRID_Y = M[t];

        TOTAL_GRID_SIZE = GRID_X * GRID_Y;
        vector<vector<map<int,int>>> grid_route(TOTAL_GRID_SIZE, vector<map<int,int>>(dFreedom));
        int count_per_curr[TOTAL_GRID_SIZE];

        gridFlat = new int[TOTAL_GRID_SIZE];

        for(int j = 0; j < TOTAL_GRID_SIZE; j++){
            gridFlat[j] = grid[t][j];
        }

        for(int j = 0; j < NODE_SIZE; j++){
            pos_i[j] = positions[t][j] / GRID_Y;
            pos_j[j] = positions[t][j] % GRID_Y;
        }

        vector<edge_route> edge(EDGE_SIZE);
        vector<vector<int>> path(NODE_SIZE);

        for (int i = 0; i < EDGE_SIZE; ++i) {
            a = vector_edges[i].first;
            b = vector_edges[i].second;
            edge[i] = edge_route(a, b, edges_cost[t][make_pair(a,b)]);
        }

        sort(edge.begin(), edge.end(), compare_edge);
        
        // print_grid_route(grid_route, GRID_X, GRID_Y, dFreedom);
        
        for (int i = 0; i < EDGE_SIZE; ++i) {
            a = edge[i].a; //vector_edges[t*EDGE_SIZE+i].first;
            b = edge[i].b; //vector_edges[t*EDGE_SIZE+i].second;
            pos_a_x = pos_i[a];
            pos_a_y = pos_j[a];
            pos_b_x = pos_i[b];
            pos_b_y = pos_j[b];

            //printf("%2d (%2d,%2d) PE: %d -> %2d (%2d,%2d) PE: %d Cost: %d\n", a, pos_a_x, pos_a_y, positions[t][a], b, pos_b_x, pos_b_y, positions[t][b], edge[i].cost);

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

#endif
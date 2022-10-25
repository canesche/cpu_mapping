#ifndef __NEW_ROUTING_H
#define __NEW_ROUTING_H

using namespace std;

#include <vector>
#include <map>
#include <algorithm>
#include <MyPathList.h>


struct edge_route {
    int a, b, cost;
    edge_route(int i = 0, int j = 0, int c = 0) : a(i), b(j), cost(c) {}
};


int find_quadrant(int pos_a_x, int pos_a_y, int pos_b_x, int pos_b_y)
{
    if(pos_a_x < pos_b_x)
    {
        if(pos_a_y < pos_b_y) return 0;
        else return 1;
    } else
    {
        if(pos_a_y < pos_b_y) return 3;
        else return 2;
    }
};
/*
  0
3   1
  2
*/
int find_direction(pair<int,int> direction)
{
    if(direction.first == 1) return 0;
    else if(direction.first == -1) return 2;
    else if(direction.second == 1) return 1;
    else return 3;
};


void print_path(vector<pair<int,int>> &path, int pe_init_x, int pe_init_y, int a, int N, int M) {
    vector<int> pe_path(path.size());
    pe_path.push_back(pe_init_x * M + pe_init_y);
    int pe_coord_x = pe_init_x;
    int pe_coord_y = pe_init_y;
    for (int i = 0; i< path.size(); ++i)
    {
        pe_coord_x += path[i].second;
        pe_coord_y += path[i].first;
        pe_path[i] = (pe_coord_x * M + pe_coord_y);
    }

    printf("\n");
    for (int i = 0; i < N*M; ++i)
    {
        if(i % M == 0) printf("\n");
        if(find(pe_path.begin(), pe_path.end(), i) != pe_path.end()) printf("x");
        else printf("#");
    }
    printf("\n");
};

bool simulate_move(vector<vector<pair<int,int>>> &paths, vector<vector<map<int,int>>> &grid, 
    int pe_init_x, int pe_init_y, int pe_end_x, int pe_end_y,
    int a, int GRID_X, int GRID_Y)
{
    vector<int> best_route(paths.size(),0);
    vector<pair<int,int>> path;
    int step;


    for(int i=0;i<paths.size();i++)
    {
        path = paths[i]; 
        step = 0;
        int pe_curr = pe_init_x * GRID_Y + pe_init_y;
        int pe_coord_x = pe_init_x;
        int pe_coord_y = pe_init_y;
        for(int j=0;j<path.size();j++)
        {
            ++step;
            int direction = find_direction(path[j]);
            if(grid[pe_curr][direction][step] == a) /* se o caminho a seguir tem uma anotação com  o nó a no tempo step */
            {
                best_route[i]+=2;
                pe_coord_x += path[j].second;
                pe_coord_y += path[j].first;
                pe_curr = (pe_coord_x) * GRID_Y + (pe_coord_y);
            }
            else if(grid[pe_curr][direction][step] == 0) /* se o caminho está livre */
            {
                best_route[i]+=1;
                pe_coord_x += path[j].second;
                pe_coord_y += path[j].first;
                pe_curr = (pe_coord_x) * GRID_Y + (pe_coord_y);
            }
            else /* se o caminho já está usado */
            {
                best_route[i] = -1;
                break;
            }
        }
    }

    int max_path_index = max_element(best_route.begin(),best_route.end()) - best_route.begin();


    if(best_route[max_path_index] == -1) 
    {
        
        return false;
    }
    else
    {
        // print_path(paths[max_path_index],pe_init_x,pe_init_y,a,GRID_X,GRID_Y);
        step = 0;
        int pe_curr = pe_init_x * GRID_Y + pe_init_y;
        int pe_coord_x = pe_init_x;
        int pe_coord_y = pe_init_y;
        for(int i=0; i<paths[max_path_index].size();i++) /* fazendo anotações definitivas */
        {

            ++step;
            int direction = find_direction(paths[max_path_index][i]);

            // if (a==26) cout << pe_curr << "->";

            grid[pe_curr][direction][step] = a; /* fazendo anotações */

            pe_coord_x += paths[max_path_index][i].second;
            pe_coord_y += paths[max_path_index][i].first;
            pe_curr = (pe_coord_x) * GRID_Y + (pe_coord_y);
        }

        // if(a==26) cout << "(" << pe_curr << ":" << pe_end_x * GRID_Y + pe_end_y <<  ")\n";
        return true;
    }

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
};


bool compare_edge(edge_route &a, edge_route &b) {
    return a.cost > b.cost;
};

void routing(
    int* N, 
    int* M, 
    int EDGE_SIZE, 
    int NODE_SIZE, 
    vector<pair<int,int>> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    int times, 
    bool* successfulRoutings, 
    int **grid, 
    int **positions,
    int &bad_route, 
    int ARCH
) {
    
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

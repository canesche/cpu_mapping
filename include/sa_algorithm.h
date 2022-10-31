#ifndef __SA_ALGORITHM__H
#define __SA_ALGORITHM__H

#include <iostream>
#include <cmath>

int gridCost(int edges, int N, int M, int *h_edgeA, int *h_edgeB, int *positions, int LIMIT, int F_ID, int ARCH);
void sa_algorithm(int nodes, int edges, int N, int M, int &cost, int *grid, int *positions, int *v_i, int *v, vector<int> &A, int arch, double *randomvec, int &swaps, int LIMIT, int F_ID, unsigned int TEMP, unsigned int NGRIDS); 
void edgesCostConstructor(Graph g, vector<map<pair<int,int>,int>> &edges_cost, int** grid, int **positions, int* N, int* M, int ARCH, unsigned int NGRIDS);

void sa_main(
    string bench,
    const int NGRIDS,
    const int ARCH,
    const int LIMIT, 
    const int F_ID, 
    const int TEMP
) {

    int *results = new int[NGRIDS];

    Graph g(bench);
    const int nodes = g.num_nodes();
    const int edges = g.num_edges();

    // Adding size of grid
    const int dim = ceil(sqrt(nodes));

    int *N = new int[NGRIDS];
    int *M = new int[NGRIDS];

    for (int i = 0; i < NGRIDS; ++i) {
        N[i] = dim; //ceil(sqrt(nodes));
        M[i] = dim; //ceil(sqrt(nodes));
    }
    
    vector<pair<int,int>> edge_list = g.get_edges();
    vector<int> A;

    double *randomvec = new double[1000000];
    for(int i = 0; i < 1000000; i++){
        randomvec[i] = (double) rand() / (double)(RAND_MAX);
    }

    int *v = new int[nodes];
    int *v_i = new int[nodes];
    int *h_edgeA = new int[edges];
    int *h_edgeB = new int[edges];

    // creating edge list
    create_edge_list(v, v_i, h_edgeA, h_edgeB, edge_list, A, nodes, edges);   

    //Vetores para guardar o numero de buffers por aresta
    vector<map<pair<int,int>,int>> edges_cost(NGRIDS);
    vector<pair<int,int>> edges_ = g.get_edges();

    for(int k = 0; k < NGRIDS; k++){
        for(int i = 0; i < edges; i++){
            edges_cost[k][edges_[i]] = 0;
        }   
    }

    double time_total = 0.0;
    int cost_min = -1;

    //Aloca memoria pros grids
    int **grid = new int*[NGRIDS];
    int **positions = new int*[NGRIDS];
    int *local_swaps = new int[NGRIDS];
    
    for(int i = 0; i < NGRIDS; i++) {
        local_swaps[i] = 0;
        grid[i] = new int[N[i]*M[i]];
        positions[i] = new int[nodes];
        for(int j = 0; j < nodes; j++) 
            positions[i][j] = 0;
    }    

    // Fill the grid and positions
    fill_grid(g, grid, positions, nodes, N, M, NGRIDS);

    int cost = -1;

    auto start = high_resolution_clock::now();

    #pragma omp parallel for
    for (int i = 0; i < NGRIDS; ++i) {
        cost = gridCost(edges, N[i], M[i], h_edgeA, h_edgeB, positions[i], LIMIT, F_ID, ARCH);

        if(cost == edges) { 
            results[i] = cost;
            continue; // perfect solution
        }

        sa_algorithm(nodes, edges, N[i], M[i], cost, grid[i], positions[i], v_i, v, A, ARCH, randomvec, local_swaps[i], LIMIT, F_ID, TEMP, NGRIDS);
        results[i] = cost;
    }
    auto stop = high_resolution_clock::now();

    std::chrono::duration<double, std::milli> duration = (stop-start);

    time_total = duration.count();

    //printGrid(grid[0], N[0], M[0]);

    //construct edges_cost
    edgesCostConstructor(g, edges_cost, grid, positions, N, M, ARCH, NGRIDS); 

    bool *successfullRoutings = new bool[NGRIDS];

    vector<string> v_name = split(bench, "/");
    string name = v_name[v_name.size()-1];

    /*
    printf("before: \n");
    print_results(edges, edge_list, edges_cost, name, local_swaps, successfullRoutings, 
              time_total, 0, N, M, grid, g, ARCH, F_ID, 0, NGRIDS);
    */
    
    int bad_route = 0;

    start = high_resolution_clock::now();
    routing(N, M, edges, nodes, edge_list, edges_cost, NGRIDS, successfullRoutings, grid, positions, bad_route, ARCH);
    duration = (high_resolution_clock::now()-start);

    double time_routing_total = duration.count(); 

    /*
    printf("after: \n");
    print_results(edges, edge_list, edges_cost, name, local_swaps, successfullRoutings, 
              time_total, time_routing_total, N, M, grid, g, ARCH, F_ID, bad_route, NGRIDS);
    */

    int best_index = get_best_index_sa(NGRIDS, edge_list, edges_cost, successfullRoutings);

    output_graph_sa(edge_list, edges_cost, best_index, g.num_edges(), g, name, "sa", NGRIDS, ARCH);  
    
    delete h_edgeA;
    delete h_edgeB;
    delete v;
    delete v_i;
    delete local_swaps;
    delete successfullRoutings;

    edges_cost.clear();
    A.clear();
}

void printGraphInfo(int nodes, int edges, int *h_edgeA, int *h_edgeB, vector<int> &A, int *v, int *v_i){
    printf("NODES: %d EDGES: %d\n", nodes, edges);
    for (int i = 0; i < edges; ++i) printf("%d -> %d\n", h_edgeA[i], h_edgeB[i]);    
    printf("v: ");
    for (int i = 0; i < nodes; ++i) printf("%d ", v[i]);
    printf("\n");
    printf("v_i: ");
    for (int i = 0; i < nodes; ++i) printf("%d ", v_i[i]);
    printf("\n");
    printf("A: ");
    for (int i = 0; i < A.size(); ++i) printf("%d ", A[i]);
    printf("\n\n");
}

void printPlacementInfo(int nodes, int gridSize, int *grid, int *positions, int cost){
    printf("grid: ");
    for (int i = 0; i < gridSize; ++i) printf("%d ", grid[i]);
    printf("\n");
    printf("positions: ");
    for (int i = 0; i < nodes; ++i) printf("%d ", positions[i]);
    printf("\n");
    printf("cost: %d", cost);
    printf("\n\n");
}

/* Heuristic that compute the distance cost (increment) of two nodes 
 * that compose a edge in GRN that was placed in the arch, the follow
 * fuctions was used as possibilites:
 * 
 * 1) linear expresion
 * 2) polynomial expresion
 * 3) piecewise-defined function: 
 *      give a limit value LIMIT, 
 *      if increment > LIMIT -> increment *= number of edges in the GRN
 *      else return increment as it is
 * 4) exponencial expresion
 */
inline int f(int increment, int edges, int LIMIT, int F_ID){
    
    unsigned int new_increment = 0;

    if (F_ID == 0) {
        return increment;
    } else if (F_ID == 1) {
        new_increment = ((increment*increment));
    } else if (F_ID == 2) {
        if(increment <= LIMIT) return increment;
        else new_increment = increment * edges;
    } else if (F_ID == 3) {
        //testar se chega no limite de 2^31
        new_increment = 1 << (increment - 1);
    } else {
        printf("Cost fuction not defined!\n");
    }

    return new_increment;
}

int gridCost(int edges, int N, int M, int *h_edgeA, int *h_edgeB, int *positions, int LIMIT, int F_ID, int ARCH){
    
    unsigned int local_cost = 0;
    int ifrom, jfrom, ito, jto, k, increment, distManhattanI, distManhattanJ;

    for(k = 0; k < edges; k++){
        ifrom = positions[h_edgeA[k]] / M; 
        //cout << "x1= " << ifrom << " ";
        jfrom = positions[h_edgeA[k]] % M; 
        //cout << "y1= " << jfrom << " ";
        ito = positions[h_edgeB[k]] / M; 
        //cout << "x2= " << ito << " "; 
        jto = positions[h_edgeB[k]] % M;
        //cout << "y2= " << jto << " ";
        distManhattanI = abs(ito - ifrom);
        distManhattanJ = abs(jto - jfrom);

        if (ARCH == 0)
            increment = tablemesh[distManhattanI][distManhattanJ];
        else if (ARCH == 1)
            increment = table1hop[distManhattanI][distManhattanJ];
        else if (ARCH == 2)
            increment = (abs(jfrom-ifrom)%2 == 0) ? tablechess1hop[distManhattanI][distManhattanJ] : tablechessmesh[distManhattanI][distManhattanJ];
        else if (ARCH == 3)
            increment = tablehex[distManhattanI][distManhattanJ];
        else {
            printf("Architecture not defined!\n");
            return -1;
        }
        local_cost += f(increment, edges, LIMIT, F_ID); // insert function here using inline, we can use define to build code in compilation time.
    }
    return local_cost;
}


inline int local_grid_cost(int *localPositions, int* v_i, vector<int> &A, int i, 
                           int node, int N, int M, int edges, int LIMIT, int F_ID, int ARCH) {
    
    int ifrom, jfrom, pos, ito, jto, distManhattanI, distManhattanJ;

    ifrom = localPositions[node]/M; 
    jfrom = localPositions[node]%M; 
    pos = A[v_i[node]+i];
    ito = localPositions[pos]/M; 
    jto = localPositions[pos]%M;
    distManhattanJ = abs(jto - jfrom);
    distManhattanI = abs(ito - ifrom);

    if (ARCH == 0)
        return f(tablemesh[distManhattanI][distManhattanJ], edges, LIMIT, F_ID);
    else if (ARCH == 1)
        return f(table1hop[distManhattanI][distManhattanJ], edges, LIMIT, F_ID);
    else if (ARCH == 2)
        return f((abs(jfrom-ifrom) % 2 == 0) ? tablechess1hop[distManhattanI][distManhattanJ] : tablechessmesh[distManhattanI][distManhattanJ], edges, LIMIT, F_ID);
    else if (ARCH == 3)
        return f(tablehex[distManhattanI][distManhattanJ], edges, LIMIT, F_ID);
    else {
        printf("Architecture not defined!\n");
        return -1;
    }
}

void sa_algorithm(int nodes, 
               int edges, 
               int N, 
               int M, 
               int &cost, 
               int *grid, 
               int *positions, 
               int *v_i, 
               int *v, 
               vector<int> &A, 
               int arch,
               double *randomvec, 
               int &swaps, 
               int LIMIT, 
               int F_ID, 
               unsigned int TEMP,
               unsigned int NGRIDS
            ) {
    
    int currentCost = cost, nextCost, swapCount = 0, increment, distManhattanI, distManhattanJ, chess;
    const int sizeGrid = N*M;
    
    int *localGrid = new int[sizeGrid];
    for(int i = 0; i < sizeGrid; ++i) localGrid[i] = grid[i];
    
    int *localPositions = new int[nodes];
    for(int i = 0; i < nodes; ++i) localPositions[i] = positions[i]; 

    //random vector index
    double random, valor;
    int randomctrl = 0, node1, node2, old1, old2, ifrom, jfrom, ito, jto, pos;
    double T = TEMP;
    int random_count = 0;

    while(T >= 0.00001){
        for(int i = 0; i < sizeGrid; i++){
            for(int j = i+1; j < sizeGrid; j++){                
                //if we're looking at 2 empty spaces, skip                   
                if(localGrid[i] == INF && localGrid[j] == INF)
                    continue;

                node1 = localGrid[i];
                node2 = localGrid[j];
                nextCost = currentCost;
            
                swaps++;
                //remove cost from object edges                
                if(node1 != INF){
                    for(int i = 0; i < v[node1]; ++i){
                        nextCost -= local_grid_cost(localPositions, v_i, A, i, node1, N, M, edges, LIMIT, F_ID, arch); 
                    }
                }
                if(node2 != INF){
                    for(int i=0; i < v[node2]; ++i){
                        nextCost -= local_grid_cost(localPositions, v_i, A, i, node2, N, M, edges, LIMIT, F_ID, arch);
                    }
                }

                //swap positions
                old1 = i;
                old2 = j;
                if(node1 != INF) localPositions[node1] = old2;
                if(node2 != INF) localPositions[node2] = old1;
                localGrid[j] = node1;
                localGrid[i] = node2;

                //recalculate cost
                if(node1 != INF){
                    for(int i=0; i < v[node1]; ++i){
                        nextCost += local_grid_cost(localPositions, v_i, A, i, node1, N, M, edges, LIMIT, F_ID, arch);                         
                    }
                }
                if(node2 != INF){
                    for(int i = 0; i < v[node2]; ++i){
                        nextCost += local_grid_cost(localPositions, v_i, A, i, node2, N, M, edges, LIMIT, F_ID, arch);  
                    }
                }
                
                //parameter for annealing probability
                valor = exp(-1*(nextCost - currentCost)/T);
                //random number between 0 and 1
                random = randomvec[randomctrl];
                ++randomctrl;
                if(randomctrl == 1000000) randomctrl = 0;

                //if cost after changes is less than before or if cost is higher but we're in the annealing probanility range, return
                if(nextCost <= currentCost || random <= valor){
                    /* COLOCAR CONTADOR PARA VERIFICAR RANDOM */
                    // if(random <= valor && nextCost > currentCost) random_count++;
                    currentCost = nextCost;
                    swapCount++;
                } else { //else, undo changes and stay with previous cost
                    if(node1 != INF) localPositions[node1] = old1;
                    if(node2 != INF) localPositions[node2] = old2;
                    localGrid[j] = node2;
                    localGrid[i] = node1;
                }
            }
            if(cost == edges) break;             
            T *= 0.999;
        }   
        //aqui T*=0.999;
    }
    

    for(int i = 0; i < sizeGrid; i++) grid[i] = localGrid[i];
    for(int i = 0; i < nodes; i++) positions[i] = localPositions[i]; 
    cost = currentCost;

    // printf("VALIDACAO: %d\n",random_count);
}

void edgesCostConstructor(Graph g, 
                          vector<map<pair<int,int>,int>> &edges_cost, 
                          int** grid, 
                          int **positions, 
                          int* N, 
                          int* M, 
                          int ARCH,
                          unsigned int NGRIDS
                        ){
    int num_nodes = g.num_nodes();
    int num_edges = g.num_edges();
    int increment = 0;
    int distManhattanI, distManhattanJ, chess;
    vector<pair<int,int>> edges = g.get_edges();
    for(int i = 0; i < NGRIDS; i++){
        for(int j = 0; j < edges.size(); j++){
            pair<int,int> aux = edges[j];
            int A = aux.first; int B = aux.second;
            int ifrom = positions[i][A] / M[i]; 
            int jfrom = positions[i][A] % M[i]; 
            int ito = positions[i][B] / M[i]; 
            int jto = positions[i][B] % M[i];

            distManhattanI = abs(ito - ifrom);
            distManhattanJ = abs(jto - jfrom);

            if (ARCH == 0)
                increment =  distManhattanI + distManhattanJ; //tablemesh[distManhattanI][distManhattanJ];
            else if (ARCH == 1)
            {
                if(distManhattanI == 0 && distManhattanJ == 0) increment = 0;
                else increment = table1hop[distManhattanI][distManhattanJ];

            }
            else if (ARCH == 2)
            {
                int chess = (abs(jfrom-ifrom)%2==0) ? 1 : 0;
                increment = (chess==1) ? tablechess1hop[distManhattanI][distManhattanJ] : tablechessmesh[distManhattanI][distManhattanJ];
            }
            else if (ARCH == 3)
                increment = tablehex[distManhattanI][distManhattanJ];
            else
                printf("Architecture not defined!\n");
            
            //edges_cost[i].insert(pair<pair<int,int>,int>(aux,increment));
            edges_cost[i][aux] = increment;
            //cout << edges_cost[i][aux] << endl;
        }
    }
}

#endif
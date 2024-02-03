#ifndef GRID
#define GRID 

#include <vector>
#include <algorithm>

// Fill the grid in random positions
inline void fill_grid(Graph& g, 
                      int **grid, 
                      int **positions, 
                      const int nodes, 
                      int *N, 
                      int *M,
                      unsigned int NGRIDS
                    ) {
    
    vector<int> elems = g.get_nodes();
    vector<vector<int>> pos(NGRIDS);
    for(int i = 0; i < NGRIDS; i++){
        for(int j = 0; j < N[i]*M[i]; j++) {
            grid[i][j] = INF;
            if (j < nodes)
                pos[i].push_back(j);
            else
                pos[i].push_back(INF);
        }
    }
    
    for (int i = 0; i < NGRIDS; ++i) {
        std::random_shuffle(pos[i].begin(), pos[i].end()); // creating random position
        for (int j = 0; j < N[i]*M[i]; ++j) {
            grid[i][j] = pos[i][j];
        }
    }

    for (int t = 0; t < NGRIDS; t++){
        for (int i = 0; i < N[t]*M[t]; i++){
            for (int j = 0; j < nodes; j++){
                if (j == grid[t][i]) {
                    positions[t][j] = i;
                }
            }
        }
    }  
}

void printGrid(int *grid, int N, int M){
    for (int i = 0; i < N*M; ++i){
        if(i % M == 0) printf("\n");
        printf("%4d", grid[i]);
    }
    printf("\n\n");
}

void create_edge_list(int *v, int *v_i, int *h_edgeA, int *h_edgeB, 
                      vector<pair<int,int>> &edge_list, vector<int> &A, 
                      const int nodes, const int edges) {

    for(int i=0; i < nodes; ++i) {
        v[i] = 0;
        v_i[i] = 0;
    }

    int n1, n2;
    //Preenche a estrutura do grafo
    for(int i = 0; i < edge_list.size(); i++){
        n1 = edge_list[i].first;
        n2 = edge_list[i].second;
        h_edgeA[i] = n1;
        h_edgeB[i] = n2;
        v[n1]++;
        if(n1!=n2) v[n2]++;
    }

    for(int i=1; i < nodes; i++){
        v_i[i] = v_i[i-1] + v[i-1];
    }

    for(int i=0; i < nodes; i++){
        for(int j=0; j < edges; j++){
            if (h_edgeA[j] != h_edgeB[j]) {
                if(h_edgeA[j] == i) A.push_back(h_edgeB[j]);
                if(h_edgeB[j] == i) A.push_back(h_edgeA[j]);
            } else {
                if(h_edgeA[j] == i) A.push_back(h_edgeB[j]);
            }
        }
    }
}

#endif 
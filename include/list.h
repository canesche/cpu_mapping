#ifndef __LIST_H
#define __LIST_H

#include <Graph.h>

void dfs(Graph g, vector<int> &aux_edges, int *visited, int dad) {

    visited[dad] = 1;
    int child;
    vector <int> children = g.get_predecessors(dad);

    random_shuffle(children.begin(), children.end());

    for (int i = 0; i < children.size(); ++i) {
        child = children[i];
        aux_edges.push_back(dad);
        aux_edges.push_back(child);
        if (!visited[child])
            dfs(g, aux_edges, visited, child);
    }
}

void dfs_position_order(Graph g, vector<pair<pair<int,int>,int>> &EDGES, 
    const int NODE_SIZE, const int times) {

    vector<int> aux_edges, outputs;
    int *visited = new int[NODE_SIZE];

    for (int t = 0; t < times; ++t) {
        aux_edges.clear();
        outputs = g.get_outputs();
        random_shuffle(outputs.begin(), outputs.end());
        memset(visited, 0, sizeof(int)*NODE_SIZE);

        for (int i = 0; i < outputs.size(); ++i)
            dfs(g, aux_edges, visited, outputs[i]);

        for (int i = 0; i < aux_edges.size(); i += 2) {
            //printf("%d -> %d\n", aux_edges[i], aux_edges[i+1]);
            EDGES.push_back(make_pair(make_pair(aux_edges[i],aux_edges[i+1]),0));
        }
    }
}

void bfs_position_order(Graph g, vector<pair<pair<int,int>,int>> &EDGES, 
    const int NODE_SIZE, const int times) {

    std::queue<int> q;
    vector<int> aux_edges, children, outputs;
    int *visited = new int[NODE_SIZE];
    int dad, child;

    for (int t = 0; t < times; ++t) {
        aux_edges.clear();
        outputs = g.get_outputs();
        memset(visited, 0, sizeof(int)*NODE_SIZE);
        random_shuffle(outputs.begin(), outputs.end());
        
        for (int i = 0; i < outputs.size(); ++i)
            q.push(outputs[i]);

        while(!q.empty()) {
            dad = q.front();
            q.pop();

            visited[dad] = 1;

            children = g.get_predecessors(dad);

            random_shuffle(children.begin(), children.end());
            
            for (int i = 0; i < children.size(); ++i) {
                child = children[i];
                aux_edges.push_back(dad);
                aux_edges.push_back(child);
                if (!visited[child]){    
                    q.push(child);
                    visited[child] = 1;
                }
            }
        }
        for (int i = 0; i < aux_edges.size(); i += 2) {
            EDGES.push_back(make_pair(make_pair(aux_edges[i],aux_edges[i+1]),0));
        }
    }
}

void bfs_critical_path(Graph g, vector<pair<pair<int,int>,int>> &EDGES, 
    const int NODE_SIZE, const int times, int *critical_path) {

    std::queue<int> q;
    vector<int> aux_edges, children, outputs;
    int *visited = new int[NODE_SIZE];
    int dad, child;

    for (int t = 0; t < times; ++t) {
        aux_edges.clear();
        outputs = g.get_outputs();
        memset(visited, 0, sizeof(int)*NODE_SIZE);
        random_shuffle(outputs.begin(), outputs.end());
        
        for (int i = 0; i < outputs.size(); ++i)
            q.push(outputs[i]);

        while(!q.empty()) {
            dad = q.front();
            q.pop();

            visited[dad] = 1;

            children = g.get_predecessors(dad);

            random_shuffle(children.begin(), children.end());
            
            for (int i = 0; i < children.size(); ++i) {
                child = children[i];
                aux_edges.push_back(dad);
                aux_edges.push_back(child);
                if (!visited[child]){    
                    q.push(child);
                    visited[child] = 1;
                }
            }
        }
        for (int i = 0; i < aux_edges.size(); i += 2) {
            EDGES.push_back(make_pair(make_pair(aux_edges[i],aux_edges[i+1]),0));
        }
    }
}

void create_list_borders(Graph g, const int NODE_SIZE, const int GRID_SIZE, 
    int *list_borders) {

    std::queue<pair<int,int>> q;
    vector<int> son, inputs;
    int dad, child, new_cost, cost, distance;
#if __ARCH == 0
    distance = max(GRID_SIZE/2,1);
#elif __ARCH == 1
    distance = max(GRID_SIZE/2-1,1);
#endif

#if __THRESHOlD_IO > 0
    distance = min(distance, __THRESHOlD_IO);
#endif


    for (int i = 0; i < NODE_SIZE; ++i) list_borders[i] = 0;
   
    inputs = g.get_inputs();
    for (int i = 0; i < inputs.size(); ++i)
        q.push(make_pair(inputs[i],0));
    
    while(!q.empty()) {
        dad = q.front().first;
        cost = q.front().second;
        q.pop();
        if (cost > distance) continue;
        if (list_borders[dad] == 0 && cost > list_borders[dad]) 
            list_borders[dad] = cost;
        else if (cost < list_borders[dad])
            list_borders[dad] = cost;

        son = g.get_sucessors(dad);
        for (int i = 0, n = son.size(); i < n; ++i) {
            child = son[i];
            if (dad == child) continue;
            q.push(make_pair(child,cost+1));
        }
    }

    inputs = g.get_outputs();
    for (int i = 0; i < inputs.size(); ++i)
        q.push(make_pair(inputs[i],0));
    
    while(!q.empty()) {
        dad = q.front().first;
        cost = q.front().second;
        q.pop();
        if (cost > distance) continue;
        if (list_borders[dad] == 0 && cost > list_borders[dad]) 
            list_borders[dad] = cost;
        else if (cost < list_borders[dad])
            list_borders[dad] = cost;

        son = g.get_predecessors(dad);
        for (int i = 0, n = son.size(); i < n; ++i) {
            child = son[i];
            if (dad == child) continue;
            q.push(make_pair(child,cost+1));
        }
    }
}

pair<int, int> func_key(int v1, int v2) {
    return make_pair(v1, v2);
}

void remove_element(vector<int> &elem, int value) {
    elem.erase(remove(elem.begin(), elem.end(), value), elem.end());
}

void create_list_zigzag(
    Graph g, 
    vector<tuple3> &EDGES, 
    vector<vector<pair<int,int>>> &CYCLE, 
    int change, 
    const int times
){
    
    vector<int> outputList;
    stack<pair<int,int>> s;
    vector<vector<int>> L_fanin, L_fanout;
    vector<double> bc;
    const int NODE_SIZE = g.get_nodes().size();
    bool* VISITED = new bool[g.num_nodes()];
    int* critical_path = new int[g.num_nodes()];
    int a, b, direction, fanin, fanout, n_random;
    vector<pair<int,int>> aux_CYCLE;

    bc = g.get_betweenness_centrality();
    get_critical_path(g, NODE_SIZE, critical_path);

    for (int t = 0; t < times; ++t) {
        aux_CYCLE.resize(0);
        outputList = g.get_outputs();
        random_shuffle(outputList.begin(), outputList.end());
        
        n_random = 0;
        if (change == 3)
            n_random = rand() % (change + 1);
        //printf("random= %d\n", n_random);

        for (int i = 0; i < g.num_nodes(); ++i) VISITED[i] = false;

        // insert ouput in top list
        for (int i=0, n=outputList.size(); i < n; ++i)
            s.push(make_pair(outputList[i],1));
        
        L_fanin = g.get_fanin();
        L_fanout = g.get_fanout();

        while (!s.empty()) {
            //a, direction = Stack.pop(0) # get the top1
            a = s.top().first;
            direction = s.top().second;
            s.pop();

            fanin = L_fanin[a].size();
            fanout = L_fanout[a].size();

            if (direction == 1) { // direction IN

                if (fanout >= 1) {
                    
                    if (n_random == 0 && fanout > 1) { // greedy neihborhood  
                        b = L_fanout[a][0];
                        for (int i = 1; i < fanout; ++i) {
                            int aux = L_fanout[a][i];
                            if (L_fanout[aux].size() > L_fanout[b].size()) {
                                b = aux;
                            }
                        }
                    } else if (n_random == 1 && fanout > 1) { // betweens centrality 
                        b = L_fanout[a][0];
                        for (int i = 1; i < fanout; ++i) {
                            if (bc[L_fanout[a][i]] > bc[b]) b = L_fanout[a][i];
                        }
                    } else if (n_random == 2 && fanout > 1) { // critical path
                        b = L_fanout[a][0];
                        for (int i = 1; i < fanout; ++i) {
                            int aux = L_fanout[a][i];
                            if (critical_path[aux] > critical_path[b]) b = aux;
                        }
                    } else if (n_random == 3 && fanout > 1) { // zigzag
                        b = L_fanout[a].back();
                    } else { // random
                        b = L_fanout[a][rand() % L_fanout[a].size()];
                    }
                    
                    for (int i = 0; i < fanin; ++i)
                        s.push(make_pair(a,1));
                    s.push(make_pair(b,0));
                    
                    remove_element(L_fanout[a], b);
                    remove_element(L_fanin[b], a);

                    
                    if (VISITED[b]) {
                        aux_CYCLE.push_back(make_pair(a,b));
                    }
                    EDGES.push_back(tuple3(a,b,0));     

                } else if (fanin >= 1) {

                    if (n_random == 0 && fanin > 1) { // greedy neihborhood  
                        b = L_fanin[a][0];
                        for (int i = 1; i < fanin; ++i) {
                            int aux = L_fanin[a][i];
                            if (L_fanin[aux].size() > L_fanin[b].size()) b = aux;
                        }
                    } else if (n_random == 1 && fanin > 1){ // betwenees centrality
                        b = L_fanin[a][0];
                        for (int i = 1; i < fanin; ++i)
                            if (bc[L_fanin[a][i]] > bc[b]) b = L_fanin[a][i]; 
                    } else if (n_random == 2 && fanin > 1) { // critical path
                        b = L_fanin[a][0];
                        for (int i = 1; i < fanin; ++i) {
                            int aux = L_fanin[a][i];
                            if (critical_path[aux] > critical_path[b]) b = aux;
                        }
                    } else if (n_random == 3 && fanin > 1) { // zigzag
                        b = L_fanin[a].back();
                    } else { // random
                        b = L_fanin[a][rand() % L_fanin[a].size()];
                    }

                    s.push(make_pair(a,1));
                    for (int i = 0; i < fanin; ++i)
                        s.push(make_pair(b,1));

                    remove_element(L_fanin[a], b); 
                    remove_element(L_fanout[b], a);
                    
                    
                    if (VISITED[b]) {
                        aux_CYCLE.push_back(make_pair(a,b));
                    }
                    EDGES.push_back(tuple3(a,b,1));
                }

            } else { // direction OUT
                if (fanin >= 1){

                    if (n_random == 0 && fanin > 1) { // greedy neihborhood     
                        b = L_fanin[a][0];
                        for (int i = 1; i < fanin; ++i) {
                            int aux = L_fanin[a][i];
                            if (L_fanin[aux].size() > L_fanin[b].size()) {
                                b = aux;
                            }
                        }
                    } else if (n_random == 1 && fanin > 1){ // betwenees centrality
                        b = L_fanin[a][0];
                        for (int i = 1; i < fanin; ++i)
                            if (bc[L_fanin[a][i]] > bc[b]) b = L_fanin[a][i];  
                    } else if (n_random == 2 && fanin > 1){ // critical path
                        b = L_fanin[a][0];
                        for (int i = 1; i < fanin; ++i) {
                            int aux = L_fanin[a][i];
                            if (critical_path[aux] > critical_path[b]) b = aux;
                        }
                    } else if (n_random == 3 && fanin > 1) { // zigzag
                        b = L_fanin[a][0];
                    } else { // random
                        b = L_fanin[a][rand() % L_fanin[a].size()];
                    }

                    for (int i = 0; i < fanout; ++i)
                        s.push(make_pair(a,0));
                    s.push(make_pair(b,1));

                    remove_element(L_fanin[a], b); 
                    remove_element(L_fanout[b], a);
                    
                    
                    if (VISITED[b]) {
                        aux_CYCLE.push_back(make_pair(a,b));
                    }
                    EDGES.push_back(tuple3(a,b,1));    

                } else if (fanout >= 1) {

                    if (n_random == 0 && fanout > 1) { // greedy neihborhood  
                        b = L_fanout[a][0];
                        for (int i = 1; i < fanout; ++i) {
                            int aux = L_fanout[a][i];
                            if (L_fanout[aux].size() > L_fanout[b].size()) b = aux;
                        }
                    } else if (n_random == 1 && fanout > 1){ // betwenees centrality
                        b = L_fanout[a][0];
                        for (int i = 1; i < fanout; ++i)
                            if (bc[L_fanout[a][i]] > bc[b]) b = L_fanout[a][i];
                    } else if (n_random == 2 && fanout > 1) { // critical path
                        b = L_fanout[a][0];
                        for (int i = 1; i < fanout; ++i) {
                            int aux = L_fanout[a][i];
                            if (critical_path[aux] > critical_path[b]) b = aux;
                        }
                    } else if (n_random == 3 && fanout > 1) { // zigzag
                        b = L_fanout[a][0];
                    } else { // random
                        b = L_fanout[a][rand() % L_fanout[a].size()];
                    }

                    s.push(make_pair(a,0));
                    for (int i = 0; i < fanout; ++i)
                        s.push(make_pair(b,0));

                    remove_element(L_fanout[a], b); 
                    remove_element(L_fanin[b], a);
                    
                    
                    if (VISITED[b]) {
                        aux_CYCLE.push_back(make_pair(a,b));
                    }
                    EDGES.push_back(tuple3(a,b,0));
                }
            }
            VISITED[a] = true;
        }
        /*
        for (int i = 0; i < aux_CYCLE.size(); ++i) {
            printf("%d,%d ", aux_CYCLE[i].first, aux_CYCLE[i].second);
        }
        printf("\n");*/
        CYCLE.push_back(aux_CYCLE);
    } 
}

void smart_transversal_algorithm(
    Graph g, 
    vector<pair<int,int>> edges, 
    const int NODE_SIZE, 
    vector<tuple3> &EDGES, 
    vector<map<pair<int, int>,vector<tuple3>>> &dic_CYCLE, 
    int change, 
    const int times
){

    vector<vector<pair<int,int>>> CYCLE;
    map<pair<int, int>,vector<tuple3>> aux_CYCLE;

    // Create the list zigzag
    create_list_zigzag(g, EDGES, CYCLE, change, times);

    /*
    for (int i = 0; i < EDGES.size(); ++i) {
        if (i % g.num_edges() == 0) printf("\n");
        printf("%d-%d, ", EDGES[i].v0, EDGES[i].v1);
    }
    printf("\n");*/
#if __CYCLE == 1

    printf("cycle.....\n");

    const int EDGE_SIZE = g.num_edges();

    pair<int, int> key;
    bool found_start;
    int count, elem_cycle_begin, elem_cycle_end, value1, value2;
    vector<pair<int, int>> walk_key;
    vector<tuple3> dic_actual;
    int elem = 0;

    for (int t = 0; t < times; ++t) {
        aux_CYCLE.clear();
        walk_key.clear();
        dic_actual.clear();
        
        for (int i = 0; i < EDGE_SIZE; ++i){
            elem = i+t*EDGE_SIZE;
            key = func_key(EDGES[elem].v0, EDGES[elem].v1);
            aux_CYCLE[key];
            //printf("%d,%d\n", key.first, key.second);
        }

        for (int i = 0, n = CYCLE[t].size(); i < n; ++i) {
            found_start = false;
            count = 0;
            elem_cycle_begin = CYCLE[t][i].first;
            elem_cycle_end = CYCLE[t][i].second;
            //printf("begin = %d, end = %d\n", elem_cycle_begin, elem_cycle_end);
            value1 = -1;
            value2 = -1;
            
            walk_key.clear();
            for (int j = EDGE_SIZE-1; j >= 0; --j) {
                elem = j+t*EDGE_SIZE;

                if (elem_cycle_begin == EDGES[elem].v1 and !found_start) {
                    value1 = EDGES[elem].v0;
                    value2 = EDGES[elem].v1;
                    //cout << value1 << " " << value2 << endl;
                    key = func_key(value1, value2);
                    aux_CYCLE[key].push_back(tuple3(elem_cycle_end, count, 1));
                    count += 1;
                    found_start = true;
                } else if (found_start && (value1==EDGES[elem].v1 || elem_cycle_end==EDGES[elem].v0)){

                    value1 = EDGES[elem].v0;
                    value2 = EDGES[elem].v1;
                    key = func_key(value1, value2);

                    if (value1 != elem_cycle_end && value2 != elem_cycle_end) {
                        walk_key.insert(walk_key.begin()+0, key); // insert like stack
                        aux_CYCLE[key].push_back(tuple3(elem_cycle_end, count, 1));
                        count += 1;
                    } else {
                        found_start = false;

                        // Go back and update values
                        for (int k = 0; k < count/2; ++k) {
                            dic_actual = aux_CYCLE[walk_key[k]];
                            for (int l = 0, nl = dic_actual.size(); l < nl; ++l){
                                if (dic_actual[l].v0 == elem_cycle_end && (k+1) != dic_actual[l].v1) {
                                    aux_CYCLE[walk_key[k]][l].v1 = k+1;
                                    aux_CYCLE[walk_key[k]][l].v2 = 0;
                                } 
                            }
                        }
                        break; // to the next on the vector CYCLE
                    }
                }
            }
        }
        dic_CYCLE.push_back(aux_CYCLE);
    }

#if __DEBUG == 1
    printf("List Cycle\n");
    for (int t = 0; t < times; ++t) {
        printf("time: %d\n", t);
        for (int i = 0; i < EDGE_SIZE; ++i) {
            key = func_key(EDGES[i+t*EDGE_SIZE].v0, EDGES[i+t*EDGE_SIZE].v1);
            printf("%d_%d: ",EDGES[i+t*EDGE_SIZE].v0, EDGES[i+t*EDGE_SIZE].v1);
            for (int j = 0; j < dic_CYCLE[t][key].size(); ++j) {
                printf("[%d,%d,%d]", dic_CYCLE[t][key][j].v0, dic_CYCLE[t][key][j].v1,dic_CYCLE[t][key][j].v2);
            }
            printf("\n");
        }
        printf("\n\n");
    }
#endif

#endif
}

#endif
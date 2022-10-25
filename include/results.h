#ifndef RESULTS
#define RESULTS

#include <set>


string f_id[4] = {"linear","poly","limiar","exp"};

struct Results {
    string name;
    int index, menor, maior, swaps, total, N, M;
    double time_placement, time_routing;
    bool routing;
    map<int,int> histogram;

    Results() {}
    Results(string name, int idx, int me, int ma, 
                 double t_place, double t_route, 
                 bool route, map<int,int> hist, int sw, 
                 int t, int n, int m) : name(name), 
                 index(idx), menor(me), maior(ma), 
                 time_placement(t_place), time_routing(t_route),
                 routing(route), histogram(hist), swaps(sw), 
                 total(t), N(n), M(m) { }
};

bool compare(Results a, Results b) {
    if (a.maior < b.maior)
        return true;
    else if (a.total < b.total && a.maior == b.maior)
        return true;
    return false;
}

void print_graph_info(Graph &g) {
    int sn = -1;
    int min_neigh = 4;
    int neigh_by_dist = 4;
    int min_dist = 1;
    int predecessors_sum = 0;

    for (int i=0; i<g.num_nodes(); ++i)
        predecessors_sum += g.get_predecessors(i).size();

    double in_degree_mean = predecessors_sum / double(g.num_nodes());
    int n_nodes = g.num_nodes();
    int n_edges = g.num_edges();

    /* verify if a predecessor is a successor
     * and get the number of unique predecessor and successor
     */
    vector<set<int>> neighbors_by_node;
    for (int i=0; i<g.num_nodes(); ++i)
    {
        set<int> aux;
        vector<int> pred = g.get_predecessors(i);
        vector<int> succ = g.get_sucessors(i);
        for(int j=0;j<pred.size();++j)
            aux.insert(pred[j]);
        for(int j=0;j<succ.size();++j)
            aux.insert(succ[j]);
        neighbors_by_node.push_back(aux);
    }

    for (int i=0; i<g.num_nodes(); ++i)
        sn = max(int(neighbors_by_node[i].size()),sn);
    
    printf("%d\n",sn);

    while(min_neigh < sn)
    {
        neigh_by_dist += 4;
        min_neigh += neigh_by_dist;
        min_dist += 1;
    }

    printf("%d/%d\n%.3lf\n%d\n",n_nodes,n_edges,in_degree_mean,min_dist);


    
}

void print_results(int edges, vector<pair<int,int>> &edge_list,
              vector<map<pair<int,int>,int>> &edges_cost,
              string name, int *local_swaps, bool* successfullRoutings, 
              double time_total, double time_routing_total,
              int *N, int *M, int **grid, Graph &g, int ARCH, int F_ID, int bad_route, int NGRIDS) {

    vector<Results> res(NGRIDS);
    map<int,int> hist;

    int a, b, sum, menor, maior, cost;
    for (int j = 0; j < NGRIDS; ++j) {
        sum = 0;
        menor = 99999;
        maior = -1;
        hist.clear();

        for (int i = 0; i < edges; ++i) {
            a = edge_list[i].first; 
            b = edge_list[i].second; 

            cost = edges_cost[j][make_pair(a,b)];

            if (hist.count(cost) > 0)
                hist[cost]++;
            else
                hist[cost] = 1;

            menor = min(menor, cost);
            maior = max(maior, cost);
            sum += cost;
        }

        res[j] = Results(name, j, menor, maior, time_total, time_routing_total, successfullRoutings[j], hist, local_swaps[j], sum, N[j], M[j]);

        /*
        ofstream f;
        if (ARCH == 0) f.open("results/"+f_id[F_ID]+"/"+name+"/mesh/"+to_string(j)+".txt", std::ofstream::out); 
        else if (ARCH == 1) f.open("results/"+f_id[F_ID]+"/"+name+"/1-hop/"+to_string(j)+".txt", std::ofstream::out);
        else if (ARCH == 2) f.open("results/"+f_id[F_ID]+"/"+name+"/chess/"+to_string(j)+".txt", std::ofstream::out);
        else if (ARCH == 3) f.open("results/"+f_id[F_ID]+"/"+name+"/hex/"+to_string(j)+".txt", std::ofstream::out);*/
        
        /*
        f << N[j] << " " << M[j] << "\n";
        printf("%d %d\n", dim, dim);
        for (int i = 0; i < N[j]*M[j]; ++i) { 
            if (grid[j][i] != INF) {
                //f << i << " " << g.get_name_node(grid[j][i]) << "\n";
                printf("%d %s\n", i, g.get_name_node(grid[j][i]).c_str());
            }
        }
        f.close();
        */
    }

    // sorting the results
    sort(res.begin(), res.end(), compare);

    int best_case = res[0].maior;
    int best_case_count = 0;
    bool just_get = false;
    for (int i=0;i<(NGRIDS-bad_route); ++i) { //get the number of solutions with the best wc
        if (res[i].maior == best_case)
            best_case_count++;
    }

    printf("%d(%d)\n", best_case, best_case_count);
    printf("(%dx%d)\n",res[0].N, res[0].M);
}

#endif
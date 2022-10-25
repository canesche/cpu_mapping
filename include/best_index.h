#ifndef BEST_INDEX
#define BEST_INDEX
 
const unsigned int get_best_index_sa (
    const int NGRIDS,
    vector<pair<int,int>> &edge_list,
    vector<map<pair<int,int>,int>> &edges_cost,
    bool *successfullRoutings
) {

    int best_index = 0;
    int sum_global = INT32_MAX;
    for (int i = 0; i < NGRIDS; ++i) {
        int sum_local = 0;
        for (int j = 0; j < edge_list.size(); ++j) {
            int a = edge_list[j].first;
            int b = edge_list[j].second;
            sum_local += edges_cost[i][make_pair(a,b)];
        }
        cout << i << " " << sum_local << endl;
        if (sum_local < sum_global && successfullRoutings[i]) {
            sum_global = sum_local;
            best_index = i;
        }
    } 
    return best_index;
}  

const unsigned int get_best_index_yott (
    const int NGRIDS,
    const int EDGE_SIZE,
    vector<tuple3> &vector_edges,
    vector<map<pair<int,int>,int>> &edges_cost,
    bool *successfullRoutings
) {

    int best_index = 0;
    int sum_global = INT32_MAX;
    int a, b, dir;
    for (int i = 0; i < NGRIDS; ++i) {
        int sum_local = 0;
        for (int j = 0; j < EDGE_SIZE; ++j) {
            a = vector_edges[i*EDGE_SIZE+j].v0;
            b = vector_edges[i*EDGE_SIZE+j].v1;
            dir = vector_edges[i*EDGE_SIZE+j].v2;
            sum_local += edges_cost[i][make_pair(a,b)];
        }
        cout << i << " " << sum_local << endl;
        if (sum_local < sum_global && successfullRoutings[i]) {
            sum_global = sum_local;
            best_index = i;
        }
    } 
    return best_index;
} 

#endif
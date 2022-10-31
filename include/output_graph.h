#ifndef __OUTPUT_GRAPH_H
#define __OUTPUT_GRAPH_H

#include <stdio.h>
#include <vector>
#include <map>
#include <Graph.h>

using namespace std;

void output_graph_yott(
    vector<tuple3> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    const int idx, 
    const int EDGE_SIZE,
    Graph &g,
    string name,
    string tool,
    int size,
    const int arch
){
    string str_arch = "";
    if (arch == 0) str_arch = "mesh";
    else if (arch == 1) str_arch = "1hop";
    else if (arch == 2) str_arch = "chess";
    else if (arch == 3) str_arch = "hex";

    int port[g.num_nodes()];
    for (int i = 0; i < g.num_nodes(); ++i) port[i] = 0;

    int a, b, dir, value;

    string output = "results/"+tool+"/"+str_arch+"/"+to_string(size)+"/"+name+".dot";

    ofstream f;
    f.open(output);
    f << "digraph g {\n";

    for (int j = 0; j < g.num_nodes(); ++j) {
        f << j << " [label = " << g.get_name_node(j).c_str() << ";op = " << g.get_name_op(j).c_str() << ";";
        //printf("%d [label = %s;op = %s;", j, g.get_name_node(j).c_str(), g.get_name_op(j).c_str());
        if (!g.get_name_value(j).empty()) {
            //printf("value = %s;",g.get_name_value(j).c_str());
            f << "value = " << g.get_name_value(j).c_str() << ";";
        }
        //printf("]\n");
        f << "]\n";
    }

    for (int j = 0; j < EDGE_SIZE; ++j) {
        a = vector_edges[idx*EDGE_SIZE+j].v0;
        b = vector_edges[idx*EDGE_SIZE+j].v1;
        dir = vector_edges[idx*EDGE_SIZE+j].v2;
        value = edges_cost[idx][make_pair(a,b)]-1;
        if (dir == 0) { 
            //printf("%d -> %d [port=%d; weight=%d;]\n", a, b, port[b]++, value);
            f << a << " -> " << b << " [port=" << port[b]++ << "; weight=" << value << ";]\n";
        } else {
            //printf("%d -> %d [port=%d; weight=%d;]\n", b, a, port[a]++, value);
            f << b << " -> " << a <<  " [port=" << port[a]++ << "; weight=" << value << ";]\n";
        } 
    }
    //printf("}\n");
    f << "}\n";
} 

void output_graph_sa(
    vector<pair<int,int>> &vector_edges, 
    vector<map<pair<int,int>,int>> &edges_cost, 
    const int idx, 
    const int EDGE_SIZE, 
    Graph &g,
    string name,
    string tool,
    int size,
    const int arch
){

    string str_arch = "";
    if (arch == 0) str_arch = "mesh";
    else if (arch == 1) str_arch = "1hop";
    else if (arch == 2) str_arch = "chess";
    else if (arch == 3) str_arch = "hex";

    int port[g.num_nodes()];
    for (int i = 0; i < g.num_nodes(); ++i) port[i] = 0;
    
    int a, b, dir, value;

    string output = "results/"+tool+"/"+str_arch+"/"+to_string(size)+"/"+name+".dot";

    ofstream f;
    f.open(output);
    f << "digraph g {\n";

    for (int j = 0; j < g.num_nodes(); ++j) {
        //printf("%d [label = %s;op = %s;", j, g.get_name_node(j).c_str(), g.get_name_op(j).c_str());
        f << j << " [label = " << g.get_name_node(j).c_str() << ";op = " << g.get_name_op(j).c_str() << ";";
        if (!g.get_name_value(j).empty()) {
            //printf("value = %s;",g.get_name_value(j).c_str());
            f << "value = " << g.get_name_value(j).c_str() << ";";
        }
        //printf("]\n");
        f << "]\n";
    }

    for (int j = 0; j < EDGE_SIZE; ++j) {
        a = vector_edges[j].first;
        b = vector_edges[j].second;
        value = edges_cost[idx][make_pair(a,b)]-1;
        //printf("%d -> %d [port=%d; weight=%d];\n", a, b, port[b]++, value);
        f << a << " -> " << b << " [port=" << port[b]++ << "; weight=" << value << ";]\n";
    }
    //printf("}\n");
    f << "}\n";
}

#endif
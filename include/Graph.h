#ifndef __GRAPH__H
#define __GRAPH__H

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>

using namespace boost;
using namespace std;

// Vertex properties
typedef boost::property <boost::vertex_name_t, string, boost::property <boost::vertex_color_t, float>> vertex_p;  
// Edge properties
typedef boost::property <boost::edge_weight_t, double> edge_p;
// Graph properties
typedef boost::property <boost::graph_name_t, string> graph_p;
// adjacency_list-based type
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS, vertex_p, edge_p, graph_p> graph_t;
// Vertex index map
typedef boost::property_map<graph_t, boost::vertex_index_t>::type VertexIndexMap;

using vertex_t = boost::graph_traits<graph_t>::vertex_descriptor;

typedef typename boost::graph_traits<graph_t>::edge_iterator   e_iter;
typedef typename boost::graph_traits<graph_t>::vertex_iterator v_iter;

class Graph {
    public:
        struct Vertex { int foo; };
        Graph();
        Graph(string filename);
        Graph(const Graph &g);
        ~Graph();
        void print();
        void print_graph_number();
        void write(string filename="test.dot");
        void collect_input(string filename, string str_input, map<int,string> &new_map);
        //vertex_t add_node(Vertex u); // ps da vida kkk
        const int num_nodes();
        void add_edge(vertex_t u, vertex_t v);
        const int num_edges();
        vector<pair<int,int>> get_edges();
        vector<pair<int,int>> get_edges_inverse();
        vector<int> get_nodes();
        string get_name_node(int u);
        string get_name_op(int u);
        string get_name_value(int u);
        vector<int> get_predecessors(int u);
        vector<vector<int>> get_fanin();
        vector<vector<int>> get_fanout();
        vector<int> get_sucessors(int u);
        vector<int> get_inputs();
        vector<int> get_outputs();
        vector<double> get_betweenness_centrality();
    private:
        graph_t graph;
        boost::dynamic_properties dp;
        vector<int> nodes;
        vector<pair<int,int>> edges;
        map<int,vector<int>> node_in_degree;
        map<int,vector<int>> node_out_degree;
        map<int,string> name_label;
        map<int,string> op_label;
        map<int,string> value_label;
};

#endif
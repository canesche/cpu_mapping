#include <Graph.h>

Graph::Graph() {
	this->graph = graph_t(0);
	this->dp = boost::dynamic_properties(boost::ignore_other_properties);
	
	boost::property_map<graph_t, boost::vertex_name_t>::type name = get(boost::vertex_name, this->graph);
	this->dp.property("node_id", name);

	boost::property_map<graph_t, boost::vertex_name_t>::type op = get(boost::vertex_name, this->graph);
	this->dp.property("op", op);

	boost::property_map<graph_t, boost::vertex_name_t>::type label = get(boost::vertex_name, this->graph);
	this->dp.property("label", label);

	boost::property_map<graph_t, boost::vertex_name_t>::type value = get(boost::vertex_name, this->graph);
	this->dp.property("value", value);

	boost::property_map<graph_t, boost::edge_weight_t>::type weight = get(boost::edge_weight, this->graph);
	this->dp.property("weight", weight);

	// Use ref_property_map to turn a graph property into a property map
	boost::ref_property_map<graph_t*, string> gname(get_property(this->graph, boost::graph_name));
	
	istringstream gvgraph("digraph { }");

	read_graphviz(gvgraph, this->graph, this->dp, "node_id");
}

Graph::Graph(string filename) {

    this->graph = graph_t(0);
	this->dp = boost::dynamic_properties(boost::ignore_other_properties);

	boost::property_map<graph_t, boost::vertex_name_t>::type name = get(boost::vertex_name, this->graph);
	this->dp.property("node_id", name);

	boost::property_map<graph_t, boost::vertex_name_t>::type label = get(boost::vertex_name, this->graph);
	this->dp.property("label", label);

	boost::property_map<graph_t, boost::vertex_name_t>::type value = get(boost::vertex_name, this->graph);
	this->dp.property("value", value);

	boost::property_map<graph_t, boost::vertex_name_t>::type op = get(boost::vertex_name, this->graph);
	this->dp.property("op", op);

	boost::property_map<graph_t, boost::edge_weight_t>::type weight = get(boost::edge_weight, this->graph);
	this->dp.property("weight", weight);

	// Use ref_property_map to turn a graph property into a property map
	//ref_property_map<graph_t*, string> gname(get_property(this->graph, graph_name));
	
	//this->dp.property("name", gname);

	if (filename.substr(filename.find_last_of(".") + 1) == "dot") {
    	ifstream dot_file = ifstream(filename);
		read_graphviz(dot_file, this->graph, this->dp, "node_id");
	} else {
		cout << "ERROR: This file " << filename << " is not a dot file!" << endl;
	}

	// set edges
	int u, v, size;
	pair<e_iter,e_iter> it;
	vector<pair<int,int>> aux_inv;
	for(it = boost::edges(this->graph); it.first != it.second; ++it.first) {
		u = source(*it.first, this->graph);
		v = target(*it.first, this->graph);
		this->edges.push_back(make_pair(u, v));
		if (u != v) node_in_degree[v].push_back(u);
		if (u != v) node_out_degree[u].push_back(v);
	}

	// set nodes
	v_iter vi, vi_end, next;
	boost::tie(vi, vi_end) = vertices(this->graph);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		name_label[*vi] = name[*vi];
		op_label[*vi] = op[*vi];
		value_label[*vi] = value[*vi];
		this->nodes.push_back(*vi);
	}
}

Graph::Graph(const Graph &g) {
	this->graph = g.graph;
	this->dp = g.dp;
	this->nodes = g.nodes;
	this->edges = g.edges;
	this->node_in_degree = g.node_in_degree;
	this->node_out_degree = g.node_out_degree;
	this->name_label = g.name_label;
	this->op_label = g.op_label;
	this->value_label = g.value_label;
}

Graph::~Graph() {
	nodes.clear();
	edges.clear();
}

void Graph::print() {
    write_graphviz_dp(std::cout, this->graph, this->dp, "node_id");
}

void Graph::write(string path) {
	string graphName;  
	if (path.length() < 4 && path.substr(path.find_last_of(".") + 1) != "dot") {
		graphName = path;
    	path.append(".dot");
	} else {
		graphName = path.substr(path.find_last_of(".") - 1);
	}

	ofstream dotfile (path.c_str());
	
	write_graphviz_dp(dotfile, this->graph, this->dp);
}

/*
vertex_t Graph::add_node(Vertex u) {
	boost::add_vertex(u, this->graph);
}
*/

const int Graph::num_nodes() {
	return this->nodes.size();
}

void Graph::add_edge(vertex_t u, vertex_t v) {
	boost::add_edge(u, v, this->graph);
}

const int Graph::num_edges() {
	return this->edges.size();
}

vector<pair<int,int>> Graph::get_edges() {
	return this->edges;
}

vector<int> Graph::get_nodes() {
	return this->nodes;
}

string Graph::get_name_node(int u) {
	return this->name_label[u];
}

string Graph::get_name_op(int u) {
	return this->op_label[u];
}

string Graph::get_name_value(int u) {
	return this->value_label[u];
}

vector<int> Graph::get_predecessors(int u) {
	return this->node_in_degree[u];
}

vector<int> Graph::get_sucessors(int u) {
	return this->node_out_degree[u];
}

vector<vector<int>> Graph::get_fanin() {
	vector<vector<int>> aux;
	for (int i = 0, n = num_nodes(); i < n; ++i)
		aux.push_back(get_predecessors(i));
	return aux;
}

vector<vector<int>> Graph::get_fanout() {
	vector<vector<int>> aux;
	for (int i = 0, n = num_nodes(); i < n; ++i)
		aux.push_back(get_sucessors(i));
	return aux;
}

vector<int> Graph::get_inputs() {
	vector<int> aux;
	for (int i = 0; i < num_nodes(); ++i)
		if (get_predecessors(i).size() == 0) aux.push_back(i);
	return aux;
}

vector<int> Graph::get_outputs() {
	vector<int> aux;
	for (int i = 0; i < num_nodes(); ++i)
		if (get_sucessors(i).size() == 0) aux.push_back(i);
	return aux;
}

vector<pair<int,int>> Graph::get_edges_inverse() {
	vector<pair<int,int>> aux;
	int u, v;
	for(int i = num_edges()-1; i >= 0; --i){
		u = get_edges()[i].first;
		v = get_edges()[i].second;
		aux.push_back(make_pair(v,u));
	}
	return aux;
}

void Graph::print_graph_number() {
	vector<pair<int,int>> edges;

	cout << "digraph G {" << endl;
	for (int i = 0; i < num_edges(); ++i) {
		cout << this->edges[i].first << "->" << this->edges[i].second << endl;
	}
	cout << "}" << endl;

}

vector<double> Graph::get_betweenness_centrality() {
	vector<double> centrality(boost::num_vertices(this->graph), 0.0);
	VertexIndexMap v_index = get(boost::vertex_index, this->graph);
	boost::iterator_property_map<vector<double>::iterator, VertexIndexMap> vertex_property_map = make_iterator_property_map(centrality.begin(), v_index);

	boost::brandes_betweenness_centrality(this->graph, vertex_property_map);

	//double max = *max_element(centrality.begin(), centrality.end());
	//max = (max > 0.0) ? max : 1;
	
	//for (int i = 0, n = centrality.size(); i < n; ++i)
	//	centrality[i] = centrality[i] / max;

	return centrality;
}
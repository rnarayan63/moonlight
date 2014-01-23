//
// Implementing Djikstra's Algorithm
//

#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <utility>
#include <queue>
#include <algorithm>
#include <random>
#include <limits>
#include <ctime>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <set>
#include <thread>
#include <sstream>

using namespace std;



// A big number to define uneachable or infinite distance
const int INF_DIST = 100000;

class minHeap {
public:
	// destination vertex and cost from source
	typedef std::pair<int, double> element;

	// key is vertex and value is index in heap
	typedef std::map<int, int> MapType;

	minHeap();
	void insert(element);
	void modify_value_lower(element e1, double cost); // find and replace e1's cost and re-heapify
	void push_up(int); // starting at this index in the heap push element as far up as possible
	void swap(int, int); // swap two elements in the heap, appropriately adjust the dictionary
	bool empty();
	element top();
	void pop();
	int find_index(int);
	int size();
	element get_item(int);
	void heapify(int);
	void print_data_index_map();
	void print_heap();
	void test_minHeap();
	~minHeap();

private:
	std::vector<element> heap;
	MapType data_index_map;
	int leftChild(int);
	int rightChild(int);
	int parent(int);

};

minHeap::minHeap() {
	// initialize heap with dummy element in index 0
	element el;
	el.first = -1;
	el.second = 0;
	heap.push_back(el);
	//std::cout << " CALLED CONSTRUCTOR " << heap.size() << std::endl;

}

minHeap::~minHeap() {
	//std::cout << " DESTRUCTOR minHeap Called" << std::endl;
	heap.resize(0);
}
void minHeap::heapify(int index) {
	using namespace std;
	static int count = 0;
	count++;
	//std::cout << "calling Heapify - " << count << " index = " << index << std::endl;
	if (heap.size() == 1)
		return;

	int l = leftChild(index);
	int r = rightChild(index);
	int min_index = index;

	// first is the destination vertex and second is the cost
	if ((l != -1) && (heap[l].second < heap[index].second))
		min_index = l;

	if ((r != -1) && (heap[r].second < heap[min_index].second))
		min_index = r;

	if (min_index != index) {
		swap(min_index, index);
		heapify(min_index);
	}
}

minHeap::element minHeap::get_item(int index) {
	return heap[index];
}

int minHeap::leftChild(int index) {
	int left = index * 2;
	if (left > size())
		return -1;
	else
		return left;
}

int minHeap::rightChild(int index) {
	int right = index * 2 + 1;
	if (right > size())
		return -1;
	else
		return right;
}
void minHeap::swap(int el1_index, int el2_index) {

	element temp = heap[el1_index];
	//update the data_index_map Map appropriately
	data_index_map[heap[el1_index].first] = el2_index;
	data_index_map[heap[el2_index].first] = el1_index;

	heap[el1_index] = heap[el2_index];
	heap[el2_index] = temp;
}

int minHeap::parent(int index) {
	if (index == 1)
		return -1;
	else
		return index / 2;
}

void minHeap::push_up(int index) {
	using namespace std;

	int par = parent(index);

	while (par != -1) {
		if (heap[index].second < heap[par].second) {
			swap(index, par);
			index = parent(index);
			par = parent(index);
		} else
			return;
	}
}

int minHeap::size() {
	//element 1 is unused, heap elements start with index 1
	return heap.size() - 1;
}

bool minHeap::empty() {
	if (heap.size() > 1)
		return false;
	else
		return true;
}

void minHeap::insert(element el) {
	heap.push_back(el);
	data_index_map[el.first] = heap.size() - 1;
	push_up(heap.size() - 1);
}

void minHeap::modify_value_lower(element el1, double cost) {

	// we assume that the new cost is always smaller - so item can only be pushed up

	// first find the index of the item in the heap
	int index = data_index_map[el1.first];

	// change its cost value
	heap[index].second = cost;

	//Now fix the heap, since the new value is guaranteed to be smaller
	// we call the push up function.
	push_up(index);
}

minHeap::element minHeap::top() {
	if (heap.size() > 1)
		return heap[1];
	else
		return heap[0];
}

int minHeap::find_index(int vertex) {
	return data_index_map[vertex];
}

void minHeap::print_data_index_map() {
	for (MapType::iterator it = data_index_map.begin();
			it != data_index_map.end(); it++) {
		std::cout << "key = " << it->first;
		std::cout << " Value = " << it->second << "\n";
	}
}
void minHeap::pop() {
	using namespace std;
	int first_index = 1;
	int last_index = heap.size() - 1;
	int popped = heap[first_index].first;

	data_index_map[heap[last_index].first] = 1;
	//print_data_index_map();
	swap(first_index, last_index);
	data_index_map.erase(popped);
	heap.resize(heap.size() - 1);

	// since element was swapped into index 1 - need to re-heapify starting there
	heapify(1);

}

void minHeap::test_minHeap() {
	int vertex[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	int cost[10] = { 8, 7, 6, 5, 3, 4, 43, 23, 13, 9 };

	element el;
	std::cout << "size of myheap = " << heap.size() << std::endl;

	for (int i = 0; i < 10; i++) {
		el.first = vertex[i];
		el.second = cost[i];
		insert(el);
		std::cout << "\nInserting vertex " << i << "Heap size " << heap.size()
												<< std::endl;
		for (int i = 1; i <= (int) heap.size(); i++) {
			el = get_item(i);
			std::cout << el.first << " " << el.second << std::endl;
		}
	}
	std::cout << "\n After all insertions done" << std::endl;
	for (int i = 1; i <= (int) heap.size(); i++) {
		el = minHeap::get_item(i);
		std::cout << el.first << " " << el.second << std::endl;
	}

	std::cout << "PRINTING THE data_index_map MAP" << std::endl;
	minHeap::print_data_index_map();

	std::cout << "index of 5 " << find_index(5) << std::endl;
	minHeap::pop();
	std::cout << "\n After Popping one item" << std::endl;
	for (int i = 1; i <= (int) heap.size(); i++) {
		el = get_item(i);
		std::cout << el.first << " " << el.second << std::endl;
	}
	print_data_index_map();
	std::cout << "index of 5 " << find_index(5) << std::endl;
}

void minHeap::print_heap() {
	element el;

	for (int i = 1; i <= (int) size(); i++) {
		el = get_item(i);
		std::cout << el.first << " " << el.second << std::endl;
	}
	std::cout << "heap size = " << size() << std::endl;

}

double generate_random_lo_hi(double lo, double hi) {
	return (lo + (double) rand() / ((double) RAND_MAX / (hi - lo)));
}

typedef enum color {
	EMPTY, RED, BLUE
} color;

// a map of destination vertex and cost
typedef std::pair<int, int> vertex_coord;
typedef std::map<std::pair<int, int>, double> EdgeMapType;
typedef std::map<std::pair<int, int>, color> EdgeColorMapType;
typedef std::pair<std::pair<int, int>, double> GraphEdgeType;
typedef std::map<vertex_coord, int> vcMapType;

// a vector of EdgeMapType for src vertex
typedef std::vector<EdgeMapType> graphEdgeMap;

// overload << for print color
std::ostream & operator <<(std::ostream & out, color s) {

	switch (s) {
	case EMPTY:
		out << ".";
		break;
	case RED:
		out << "R";
		break;
	case BLUE:
		out << "B";
		break;
	default:
		out << ".";
		break;
	}
	return out;
}

// an array containg to cost to reach a destination index node from source
typedef std::vector<double> Dest_vector;

// an array which stores the index of the previous node in the path to reach destination index node
typedef std::vector<int> Prev_vector;

class Graph {
public:
	// edge contains dst and cost - the array index is the src
	class Edge {
	public:
		Edge() :
			  x(-1), y(-1), edge_color(EMPTY),dst(INF_DIST),cost(0.0) {
		}
		;
		int x, y;
		color edge_color;
		unsigned int dst;
		std::pair<vertex_coord, vertex_coord> edge_coord;
		double cost;
	};

	// graph an array of list of edges
	typedef std::vector<std::list<Edge>> graph;

	// a list of indices in order from source to dest
	typedef std::list<Edge> Path;

	Graph(int, bool);
	void test_graph();
	void add_edge(int src, int dst, double cost, enum color c = EMPTY);
	void auto_generate_edges(double edge_density, double cost_range_lo,
			double cost_range_hi);
	void print_graph();
	int get_num_edges();
	int get_num_vertices();
	int get_v_num(vertex_coord vc);
	vertex_coord get_v_coord(int v);
	void add_vertex(vertex_coord vc);
	bool does_edge_exist(int src, int dst);
	void set_num_vertices(int);
	void get_list_of_adjacent_nodes(int, Path &);
	void get_shortest_path_to_all(int src, Dest_vector&, Prev_vector &,
			enum color c = EMPTY);
	int generate_graph_from_file(std::string filename);
	~Graph();
private:
	graph g;
	unsigned int num_vertices;
	unsigned int num_edges;
	bool directed;
	vector<vertex_coord> gv;
	vcMapType g_vc_map;
	graphEdgeMap g_e_map; //  unused
};

//Constructor
Graph::Graph(int v_count = 0, bool directed = false) {
	num_vertices = v_count;
	this->directed = directed;
	num_edges = 0;
	g.resize(v_count);
	gv.resize(v_count);
	g_e_map.resize(v_count);
}

Graph::~Graph() {
	//std::cout << " Graph Destructor CALLED" << std::endl;
	for (int i = 0; i < (int)g.size(); i++)
		g[i].resize(0);
	g.resize(0);
	g_e_map.resize(0);
	gv.resize(0);
}

void Graph::add_vertex(vertex_coord vc) {
	num_vertices++;
	g.resize(num_vertices);
	g_e_map.resize(num_vertices);
	gv.push_back(vc);
	g_vc_map[vc] = num_vertices - 1;
}
int Graph::get_v_num(vertex_coord vc) {
	return (g_vc_map[vc]);
}
vertex_coord Graph::get_v_coord(int v) {
	return gv[v];
}
void Graph::set_num_vertices(int vertices) {
	this->num_vertices = vertices;
}

bool Graph::does_edge_exist(int src, int dst) {
	return false;
}

int Graph::generate_graph_from_file(std::string filename){
	using namespace std;
	ifstream in;

	in.open(filename, ios::in);
	if (!in){
		cout << "File " << filename << " not found" << endl;
		return -1;
	}

	int num_vertices,src,dst;
	//int num_edges = 0;
	double cost = 0;

	g.clear();
	g_e_map.clear();

	//  read in first line - to get number of vertices
	in >> num_vertices;

	g.resize(num_vertices);
	g_e_map.resize(num_vertices);
	set_num_vertices(num_vertices);

	//cout << "num_vertices = " << num_vertices << endl;


	while (!in.eof()){
		in >> src >> dst >> cost;
		add_edge(src,dst,cost);
	}
	in.close();
	//cout << " |V| = " << get_num_vertices() << " |E| = " << get_num_edges() << endl;
	return get_num_vertices();

}


void Graph::get_list_of_adjacent_nodes(int src, Path & edges) {
	for (auto iter = g[src].begin(); iter != g[src].end(); iter++)
	edges.push_back(*iter);
}

void Graph::get_shortest_path_to_all(int src, Dest_vector & dest_vector,
		Prev_vector & prev_vector, enum color c) {

	dest_vector.clear();
	dest_vector.resize(num_vertices, INF_DIST);

	prev_vector.clear();
	prev_vector.resize(num_vertices, -1);

	double current_cost = 0;
	double new_cost = 0;

	prev_vector[src] = src;

	minHeap::element top, tmp1;

	dest_vector[src] = 0;

	minHeap edge_heap;
	minHeap::element edge_element;

	edge_element.first = src;
	edge_element.second = 0;
	edge_heap.insert(edge_element);

	for (int i = 0; i < (int)num_vertices; i++) {
		if (i != src) {
			edge_element.first = i;
			edge_element.second = INF_DIST;
			edge_heap.insert(edge_element);
		}

	}

	while (!edge_heap.empty()) {

		// pop the top element
		top = edge_heap.top();
		edge_heap.pop();

		current_cost = top.second;

		for (std::list<Edge>::iterator i = g[top.first].begin();
				i != g[top.first].end(); i++) {
			//
			if ((c != EMPTY) && (i->edge_color == c)) {
				new_cost = current_cost + i->cost;

				if (new_cost < dest_vector[i->dst]) {
					tmp1.first = i->dst;
					tmp1.second = dest_vector[i->dst];

					edge_heap.modify_value_lower(tmp1, new_cost);

					//store the previous vertex and the cost in the prev_vector and dest_vector
					prev_vector[i->dst] = top.first;
					dest_vector[i->dst] = new_cost;
				}
			}

		} // for loop

	} // while
}

void Graph::add_edge(int src, int dst, double cost, enum color c) {
	using namespace std;
	Edge edge;
	edge.dst = dst;
	edge.cost = cost;
	edge.edge_color = c;
	std::pair<int, int> xedge;

	xedge.first = src;
	xedge.second = dst;

	if (!g_e_map[src][xedge]) {
		g[src].push_back(edge);

		g_e_map[src][xedge] = cost;
		num_edges += 1;
	}
	xedge.first = dst;
	xedge.second = src;
	if (!directed) {
		if (!g_e_map[dst][xedge]) {
			edge.dst = src;
			g[dst].push_back(edge);
			g_e_map[dst][xedge] = cost;
			num_edges += 1;
		}
	}

}
// takes an empty graph with num_vertices
// auto generates edges based on edge_density
// auto generates cost based on lo and hi
void Graph::auto_generate_edges(double edge_density, double cost_range_lo,
		double cost_range_hi) {

	double x;
	double xcost;

	srand((unsigned int) time(0));

	Edge edge;
	bool visited[num_vertices][num_vertices];

	for (unsigned int i = 0; i < num_vertices; i++)
		for (unsigned int j = 0; j < num_vertices; j++)
			visited[i][j] = false;

	for (unsigned int i = 0; i < num_vertices; i++) {
		for (unsigned int j = 0; j < num_vertices; j++) {

			if ((i == j) || visited[i][j])
				continue;

			visited[i][j] = true;
			visited[j][i] = true;

			x = rand() % 100;
			if (x < edge_density) {

				xcost = generate_random_lo_hi(cost_range_lo, cost_range_hi);
				int c = 1 + rand() % 2;
				add_edge(i, j, xcost, (enum color) c);
			}
		}
	}
}

void Graph::test_graph() {

	double density;
	unsigned int num_vertices;

	std::cout << "Enter number of vertices off the graph - ";
	std::cin >> num_vertices;

	std::cout << "Enter the density of the graph - ";
	std::cin >> density;

	auto_generate_edges(density, 1.0, 10.0);

	print_graph();

	get_num_edges();
}

// overload << for print_graph
std::ostream & operator<<(std::ostream & out, std::list<Graph::Edge> l_e) {

	for (std::list<Graph::Edge>::iterator it = l_e.begin(); it != l_e.end();
			it++) {
		out << "(" << it->dst << "," << (int) it->cost << ")" << it->edge_color
				<< " ";
	}
	return out;
}

void Graph::print_graph() {
	for (unsigned int it = 0; it < g.size(); it++) {
		std::cout << it << "->" << g[it] << std::endl;
	}
}

int Graph::get_num_vertices() {
	return num_vertices;
}

int Graph::get_num_edges() {
	return num_edges / 2;
}

class mst {
public:
	mst(Graph *g, enum color edge_color = EMPTY, int start = 0);
	void create_mst_graph_prim(Graph *g, enum color edge_color = EMPTY,
			int start = 0);
	void print_mst_edges();
	void print_mst_summary();
	void mst_get_vertex_set(set<int> & v_set);
	~mst();
private:
	enum color mst_edge_color;
	int mst_start;
	EdgeMapType mst_map;
	EdgeColorMapType mst_edge_color_map;
	double mst_cost;
	int num_vertices;
};

mst::mst(Graph *g, enum color edge_color, int start) {
	mst_cost = 0;
	mst_start = start;
	mst_edge_color = edge_color;
	mst_map.clear();
	mst_edge_color_map.clear();
	num_vertices = 0;
	create_mst_graph_prim(g, edge_color, start);

}
mst::~mst() {
	mst_map.clear();
	mst_edge_color_map.clear();
}
void mst::print_mst_summary() {
	std::cout << "\nMST num vertices = " << num_vertices
			<< "\nMST- num edges   = " << mst_map.size()
			<< "\nMST color         = " << mst_edge_color
			<< "\nMST cost         = " << mst_cost << std::endl;
}
void mst::print_mst_edges() {
	for (auto i = mst_map.begin(); i != mst_map.end(); i++) {
		//int src, dst;
		//src = (i->first).first;
		//dst = (i->first).second;
		color e_color = mst_edge_color_map[i->first];
		std::cout << (i->first).first << " " << (i->first).second << " "
				<< i->second << " " << e_color << " " << endl;
	}
}

void mst::mst_get_vertex_set(set<int> & v_set) {
	for (auto i = mst_map.begin(); i != mst_map.end(); i++) {
		int src, dst;
		src = (i->first).first;
		dst = (i->first).second;
		v_set.insert(src);
		v_set.insert(dst);
	}
}
class compareGraphEdges {
public:
	bool operator()(GraphEdgeType &a1, GraphEdgeType &a2) {
		if (a2.second < a1.second)
			return true;
		return false;
	}
};

typedef std::priority_queue<GraphEdgeType, std::vector<GraphEdgeType>,
		compareGraphEdges> edge_pq_type;

void mst::create_mst_graph_prim(Graph *g, enum color c, int start_vertex) {

	int V = g->get_num_vertices();
	std::list<Graph::Edge> edges;
	edge_pq_type pq;
	GraphEdgeType ge_tmp;
	int mst_size = 0;
	std::vector<bool> node_added(V, false);
	std::pair<int, int> xedge;

	int start = start_vertex;
	int new_node;
	node_added[start] = true;
	mst_size = 1;
	edges.clear();
	g->get_list_of_adjacent_nodes(start, edges);
	//cout << "In MST - start = " << start << " number of adjacent nodes = " << edges.size() << endl;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		if ((*iter).edge_color == c) {
			ge_tmp.first.first = start;
			ge_tmp.first.second = (*iter).dst;
			ge_tmp.second = (*iter).cost;
			pq.push(ge_tmp);
		}

	}

	while (!pq.empty() && mst_size < V) {
		// get the edge with minimum cost
		// get the top of pq
		// See if both vertices are in the tree
		ge_tmp = pq.top();
		int src, dst;
		double xcost;
		src = ge_tmp.first.first;
		dst = ge_tmp.first.second;
		xcost = ge_tmp.second;
		pq.pop();
		if (node_added[src] && node_added[dst])
			// adding this edge will form a loop - so skip
			continue;
		else {
			if (!node_added[src]) {
				node_added[src] = true;
				mst_cost += xcost;
				mst_size++;
				new_node = src;
			} else {
				node_added[dst] = true;
				mst_cost += xcost;
				mst_size++;
				new_node = dst;
			}
			// add this edge to the mst
			xedge.first = src;
			xedge.second = dst;
			mst_map[xedge] = xcost;
			mst_edge_color_map[xedge] = c;
			//  add the edges connected to this new node to the priority queue
			edges.clear();
			g->get_list_of_adjacent_nodes(new_node, edges);
			for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
				if ((*iter).edge_color == c) {
					ge_tmp.first.first = new_node;
					ge_tmp.first.second = (*iter).dst;
					ge_tmp.second = (*iter).cost;
					pq.push(ge_tmp);
				}
			}
		}
	} // While loop
	num_vertices = mst_size;
}
;

//  Class shortest Path
class Shortest_path_from_single_src {
public:
	typedef std::list<int> Path;
	Shortest_path_from_single_src(Graph *g, int src, enum color c = EMPTY);
	void get_path(int dst, Path & path);
	double get_path_cost(int dst);
	~Shortest_path_from_single_src();

private:
	Graph * g;
	int src;
	Dest_vector dest_vector;
	Prev_vector prev_vector;
};

Shortest_path_from_single_src::Shortest_path_from_single_src(Graph *g_ptr,
		int src, enum color c) {
	g = g_ptr;
	this->src = src;

	g->get_shortest_path_to_all(src, dest_vector, prev_vector, c);
}
;

Shortest_path_from_single_src::~Shortest_path_from_single_src() {
	//std::cout << "shortest Path destructor" << std::endl;
}
;

void Shortest_path_from_single_src::get_path(int dst, Path & path) {

	path.clear();
	path.push_back(dst);
	int temp = dst;
	while ((temp != src) && prev_vector[temp] != -1) {
		path.push_front(prev_vector[temp]);
		temp = prev_vector[temp];
	}
	std::cout << " \nPath from " << src << " to " << dst << " is with cost "
			<< dest_vector[dst] << std::endl;
	for (std::list<int>::iterator l_index = path.begin(); l_index != path.end();
			l_index++)
		std::cout << " " << *l_index;
	std::cout << std::endl << std::endl;
}
;

double Shortest_path_from_single_src::get_path_cost(int dst) {
	return dest_vector[dst];

}

typedef enum hex_pos_t {
	L_T_CORNER,
	R_T_CORNER,
	L_B_CORNER,
	R_B_CORNER,
	L_EDGE,
	R_EDGE,
	T_EDGE,
	B_EDGE,
	INTERIOR
} hex_pos_t;
typedef enum class hex_state {
	EMPTY, RED, BLUE
} hex_state;

std::ostream & operator <<(std::ostream & out, hex_state s) {

	switch (s) {
	case hex_state::EMPTY:
		out << ".";
		break;
	case hex_state::RED:
		out << "R";
		break;
	case hex_state::BLUE:
		out << "B";
		break;
	default:
		out << ".";
		break;
	}
	return out;
}

enum color hex_state_to_edge(hex_state hc) {
	switch (hc) {
	case hex_state::EMPTY: {
		return EMPTY;
	}
	case hex_state::BLUE: {
		return BLUE;
	}
	case hex_state::RED: {
		return RED;
	}
	default:
		return EMPTY;

	}
}
std::ostream & operator <<(std::ostream & out, hex_pos_t s) {

	switch (s) {
	case L_EDGE:
		out << "L";
		break;
	case R_EDGE:
		out << "R";
		break;
	case T_EDGE:
		out << "T";
		break;
	case B_EDGE:
		out << "B";
		break;
	case L_T_CORNER:
		out << "LT";
		break;
	case R_T_CORNER:
		out << "RT";
		break;
	case L_B_CORNER:
		out << "LB";
		break;
	case R_B_CORNER:
		out << "RB";
		break;
	case INTERIOR:
		out << "I";
		break;

	default:
		out << ".";
		break;
	}
	return out;
}

typedef vector<int> arr_type;

class v_num {
public:
	v_num(int val) :
		v(val) {
	}
	;
	int v;
};

std::ostream & operator <<(std::ostream & out, v_num v) {

	if (v.v < 10)
		out << "  " << v.v << " ";
	else if (v.v < 100)
		out << " " << v.v << " ";
	else
		out << v.v << " ";

	return out;
}

std::pair<int, int> v_num_to_x_y(int v, int x_size, int y_size) {
	vertex_coord coord;
	coord.first = v % y_size;
	coord.second = v / x_size;
	return coord;
}

int x_y_to_v_num(int x, int y, int x_size, int y_size) {
	return (y * x_size + x);
}

class Grid {
public:
	Grid(unsigned int v);
	Grid(const Grid& xgrid);
	void print_grid();
	bool set(vertex_coord vc, hex_state c);
	Graph gr;
	int get_dimension(){return dimension;}
	void get_adjacent(vertex_coord vc, list<vertex_coord> & l_vc);
	hex_pos_t get_pos_type(vertex_coord);
	hex_state get_hex_state(vertex_coord);
	int get_num_pos_empty(){return (dimension * dimension) - pos_filled;}
	vertex_coord get_random_x_y_for_grid();
private:
	unsigned int dimension;
	hex_state **g;
	hex_pos_t **g_pos;
	int pos_filled;
	double edge_cost;

};

Grid::Grid(const Grid & xgrid){
	//cout << "COPY CONSTRUTCTOR for Grid CALLED " << endl;
	if (xgrid.dimension == 0){
		//cout << "COPY CONSTRUTCTOR for Grid (dimension 0) CALLED " << endl;
		g = NULL;
		g_pos = NULL;
		pos_filled = 0;
		dimension = 0;
		edge_cost = 1.0;
	}
	else {
		//cout << "COPY CONSTRUTCTOR for Grid CALLED " << endl;
		dimension = xgrid.dimension;
		pos_filled = xgrid.pos_filled;
		edge_cost = xgrid.edge_cost;
		gr = xgrid.gr;

		g = new hex_state *[dimension];
		for (int i = 0; i < (int)dimension; i++) {
			g[i] = new hex_state[dimension];
			for (int j = 0; j < (int)dimension; j++) {
				g[i][j] = xgrid.g[i][j];
			}
		}
		g_pos = new hex_pos_t *[dimension];
		for (int i = 0; i < (int)dimension; i++) {
			g_pos[i] = new hex_pos_t[dimension];
			for (int j = 0; j < (int)dimension; j++) {
				g_pos[i][j] = xgrid.g_pos[i][j];
			}
		}
	}
}
// Constructor - create a NxN matrix of hex_state = EMPTY
Grid::Grid(unsigned int v) {
	dimension = v;
	pos_filled = 0;
	edge_cost = 1.0;
	g = new hex_state *[dimension];
	for (int i = 0; i < (int)dimension; i++) {
		g[i] = new hex_state[dimension];
		for (int j = 0; j < (int)dimension; j++) {
			g[i][j] = hex_state::EMPTY;
		}
	}

	g_pos = new hex_pos_t *[dimension];
	for (int i = 0; i < (int)dimension; i++) {
		g_pos[i] = new hex_pos_t[dimension];
		for (int j = 0; j < (int)dimension; j++) {
			g_pos[i][j] = INTERIOR;
		}
	}
	g_pos[0][0] = L_T_CORNER;
	g_pos[0][dimension - 1] = L_B_CORNER;
	g_pos[dimension - 1][0] = R_T_CORNER;
	g_pos[dimension - 1][dimension - 1] = R_B_CORNER;

	for (int j = 1; j < (int)dimension - 1; j++)
		g_pos[0][j] = L_EDGE;
	for (int j = 1; j < (int)dimension - 1; j++)
		g_pos[dimension - 1][j] = R_EDGE;
	for (int i = 1; i < (int)dimension - 1; i++)
		g_pos[i][0] = T_EDGE;
	for (int i = 1; i < (int)dimension - 1; i++)
		g_pos[i][dimension - 1] = B_EDGE;

}

hex_pos_t Grid::get_pos_type(vertex_coord vc) {
	return g_pos[vc.first][vc.second];
}
hex_state Grid::get_hex_state(vertex_coord vc) {
	return g[vc.first][vc.second];
}
bool Grid::set(vertex_coord vc, hex_state c) {
	if (g[vc.first][vc.second] == hex_state::EMPTY) {

		g[vc.first][vc.second] = c;
		gr.add_vertex(vc);
		// now figure out edges to add
		list<vertex_coord> l_vc;
		//cout << " calling get adjacent for " << vc.first << " "  << vc.second << endl;
		get_adjacent(vc, l_vc);

		//cout << " get adjacent returned list of size " << l_vc.size() << endl;
		for (auto iter = l_vc.begin(); iter != l_vc.end(); iter++) {
			int x = (*iter).first;
			int y = (*iter).second;
			//cout << g[x][y] << "  " << c << endl;
			if (g[x][y] == c) {
				//cout << "iter " << gr.get_v_num(*iter) << endl;
				//cout << "vertex " << gr.get_v_num(vc) << endl;
				gr.add_edge(gr.get_v_num(*iter), gr.get_v_num(vc), 1.0,
						hex_state_to_edge(c));
			}
		}
		pos_filled ++;
		return true;
	}
	//cout << " Exiting Grid Set with Error " << "color found = " << g[vc.first][vc.second] << endl;
	return false;
}

void Grid::get_adjacent(vertex_coord vc, list<vertex_coord> & adj_list) {

	vertex_coord tmp;
	switch (g_pos[vc.first][vc.second]) {
	case L_T_CORNER: {
		tmp.first = 1;
		tmp.second = 0;
		adj_list.push_back(tmp);
		tmp.first = 0;
		tmp.second = 1;
		adj_list.push_back(tmp);
		break;
	}
	case R_T_CORNER: {
		tmp.first = dimension - 2;
		tmp.second = 0;
		adj_list.push_back(tmp);
		tmp.first = dimension - 1;
		tmp.second = 1;
		adj_list.push_back(tmp);
		tmp.first = dimension - 2;
		tmp.second = 1;
		adj_list.push_back(tmp);
		break;
	}
	case L_B_CORNER: {
		tmp.first = 0;
		tmp.second = dimension - 2;
		adj_list.push_back(tmp);
		tmp.first = 1;
		tmp.second = dimension - 1;
		adj_list.push_back(tmp);
		tmp.first = 1;
		tmp.second = dimension - 2;
		adj_list.push_back(tmp);
		break;
	}
	case R_B_CORNER: {
		tmp.first = dimension - 2;
		tmp.second = dimension - 1;
		adj_list.push_back(tmp);
		tmp.first = dimension - 1;
		tmp.second = dimension - 2;
		adj_list.push_back(tmp);
		break;
	}
	case T_EDGE: {
		tmp.first = vc.first;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		tmp.first = vc.first - 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		tmp.first = vc.first - 1;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		break;
	}
	case B_EDGE: {
		tmp.first = vc.first - 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		break;
	}
	case L_EDGE: {
		tmp.first = vc.first;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		break;
	}
	case R_EDGE: {
		tmp.first = vc.first;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first - 1;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first - 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		break;
	}
	case INTERIOR: {
		tmp.first = vc.first - 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second;
		adj_list.push_back(tmp);
		tmp.first = vc.first;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first + 1;
		tmp.second = vc.second - 1;
		adj_list.push_back(tmp);
		tmp.first = vc.first - 1;
		tmp.second = vc.second + 1;
		adj_list.push_back(tmp);
		break;
	}
	default:
		break;
	}
}

vertex_coord Grid::get_random_x_y_for_grid(){
	int dim = get_dimension();
	static std::default_random_engine generator;
	static std::uniform_int_distribution<int> distribution(0,dim*dim-1);
	//int rand_num = rand() % (dim * dim);
	int rand_num = distribution(generator);
	vertex_coord vc;
	vc.first = rand_num / dim;
	vc.second = rand_num % dim;
	//cout << "Random number = " << rand_num << endl;
	return (vc);

}
void Grid::print_grid() {

	int line_no = 0;
	int row = 1;
	char c = 'A';
	int size = dimension;

	//Header characters
	cout << "    ";
	for (int i = 0; i < size; i++)
		cout << static_cast<char>(c + i) << "   ";
	cout << endl;

	for (int j = 0; j < size; j++) {
		// print space if needed
		for (int k = 0; k < line_no; k++)
			cout << " ";
		cout << v_num(row);
		for (int i = 0; i < size; i++) {

			//cout << g_pos[i][j];

			cout << g[i][j];
			if (i != size - 1)
				cout << " - ";
		}
		cout << endl;
		line_no++;
		row++;
		if (j < size - 1) {
			cout << "    ";
			for (int k = 0; k < line_no; k++)
				cout << " ";
			cout << "\\ ";
			for (int k = 0; k < size - 1; k++)
				cout << "/ \\ ";
			cout << endl;
			line_no++;
		}
	}
	//Footer characters

	for (int i = 0; i < 2 * size + 2; i++)
		cout << " ";
	for (int i = 0; i < size; i++)
		cout << static_cast<char>(c + i) << "   ";
	cout << endl;

}

void test_read_graph_from_file(){
	Graph *g = new Graph(0);
	if (g->generate_graph_from_file("mst_data.txt") < 0)
		std::cout << "\nFile NOT found or Incorrect format\n" << std::endl;
	Shortest_path_from_single_src *sp = new Shortest_path_from_single_src (g, 0);
	delete g;
	delete sp;
	g->print_graph();
}


void hw2(int num_vert = 0, double graphDensity = 0) {
	using namespace std;
	Dest_vector dest_vector;
	Prev_vector prev_vector;
	std::list<int> path;
	double density;
	unsigned int num_vertices;
	std::cout << std::setprecision(2) << std::fixed;
	if (num_vert == 0) {
		std::cout << "Enter number of vertices of the graph - ";
		std::cin >> num_vertices;
	} else
		num_vertices = num_vert;

	if (graphDensity == 0) {
		std::cout << "Enter the density of the graph - ";
		std::cin >> density;
	} else
		density = graphDensity;

	if (num_vertices <= 0) {
		cout << " \nError in input - number of vertices <= 0 " << endl;
		return;
	}

	Graph *g = new Graph(num_vertices);

	//int src = 0;
	//int dst = 0;

	// generate the graph
	g->auto_generate_edges(density, 1.0, 10.0);

	g->print_graph();
	// get shortest path to all nodes from this src node 0
	Shortest_path_from_single_src *sp = new Shortest_path_from_single_src(g, 0,
			BLUE);

	// calculate total cost of all edges
	double cost_dst_sum = 0;
	for (int i = 1; i <= (int)num_vertices; i++) {
		//sp->get_path(i, path);
		cost_dst_sum += sp->get_path_cost(i);
	}

	int num_edges = g->get_num_edges();
	cout << " Requested number of vertices |V|  = " << num_vertices << endl;

	cout << " Number of edges generated |E|     = " << g->get_num_edges()
											<< endl;

	double max_edges = (num_vertices * (num_vertices - 1) / 2);
	cout << " Requested density                 = " << density << endl;
	cout << " Actual Density                    = "
			<< (num_edges / max_edges) * 100 << std::endl;
	cout << "\nTotal cost of all shortest Paths   = " << cost_dst_sum << endl;
	if (cost_dst_sum != 0)
		cout << " Average Shortest path             = "
		<< cost_dst_sum / (num_vertices - 1) << endl;
	else
		cout << " OOPs!!! No Paths exist from this Node             = " << endl;

	delete g;
	delete sp;
	return;
}


void hw3() {

	// Create a graph by reading input from a file
	// generate Minimum Spanning tree using Prim Algorithm.
	using namespace std;

	Graph *g = new Graph(0, false);

	if (g->generate_graph_from_file("mst_data.txt") < 0){
		cout << "\n File not Found or incorrect format \n " << endl;
		return;
	}

	// create mst object with this graph
	mst * test_mst = new mst(g);
	test_mst->print_mst_summary();

	delete (g);
	delete(test_mst);
}

vertex_coord fill_grid(Grid &grid, hex_state start_player_color){
	//int comp_x,comp_y;
	vertex_coord vc,vc_first;
	bool first = true;
	//int size = grid.get_dimension();
	int total = grid.get_num_pos_empty();
	hex_state next_hex_state = start_player_color;

	while (total > 0){

		start_player_color = next_hex_state;
		// keep toggling colors to simulate player versus opponent
		if (start_player_color == hex_state::BLUE)
			next_hex_state = hex_state:: RED;
		else
			next_hex_state = hex_state:: BLUE;
		
		// generate random x and y coordinates.
		vc = grid.get_random_x_y_for_grid();

		while (!grid.set(vc,start_player_color)) {

			vc = grid.get_random_x_y_for_grid();

		}

		if (first)
			vc_first = vc;

		total -= 1;
	}
	return vc_first;
}

hex_state who_won(Grid & grid, bool full_grid = true){

	int size = grid.get_dimension();

	set<int> mst_red_set;
	set<int> mst_blue_set;
	vertex_coord vc_tmp;
	hex_pos_t tmp_pos;
	bool found = false;

	// walk thru all the left EDGE grid points
	// if the position is occupied by BLUE
	// using that as the starting location - construct an MST
	// if the resulting MST contains a node from the RIGHT Edge - we are done
	// othewise go to the next node on the left EDGE.

	for (int j =  0; j < size; j++) {
		vc_tmp.first = 0;
		vc_tmp.second = j;
		tmp_pos = grid.get_pos_type(vc_tmp);
		if (grid.get_hex_state(vc_tmp) == hex_state::RED){
		//cout << *iter << " " << tmp_pos << " ";
			// create mst object with this graph
			mst *test_mst_red = new mst(&grid.gr, RED, grid.gr.get_v_num(vc_tmp));
			//test_mst_red->print_mst_edges();
			//test_mst_red->print_mst_summary();
			// grab the set of vertices in the MST
			test_mst_red->mst_get_vertex_set(mst_red_set);

			for (auto iter = mst_red_set.begin(); iter != mst_red_set.end(); iter++) {

				vc_tmp = grid.gr.get_v_coord(*iter);
				tmp_pos = grid.get_pos_type(vc_tmp);

				if ((tmp_pos == R_EDGE) || (tmp_pos == R_T_CORNER)
						|| (tmp_pos == R_B_CORNER)){
					found = true;
					break;
				}
			} // for statement
			delete (test_mst_red);
			if (found) {
				//cout << "\n\nBLUE won the GAME  !!!\n\n";
				return hex_state::RED;
			}
		}

	} // outer for loop

	if (!full_grid){
		set<int> mst_blue_set;

		found = false;

		// walk thru all the Top EDGE grid points
		// if the position is occupied by BLUE
		// using that as the starting location - construct an MST
		// if the resulting MST contains a node from the Bottom Edge - we are done
		// othewise go to the next node on the Top EDGE.
		for (int i =  0; i < size; i++){
			vc_tmp.first = i;
			vc_tmp.second = 0;
			//cout << *iter << " " << tmp_pos << " ";
			if (grid.get_hex_state(vc_tmp) == hex_state::BLUE){
				// create mst object with this graph
				mst *test_mst_blue = new mst(&grid.gr,
											BLUE,
											grid.gr.get_v_num(vc_tmp));
				//test_mst_blue->print_mst_edges();
				//test_mst_blue->print_mst_summary();

				// grab the set of vertices in the MST
				test_mst_blue->mst_get_vertex_set(mst_blue_set);
				//int num_left = 0;
				//int num_right = 0;

				for (auto iter = mst_blue_set.begin(); iter != mst_blue_set.end(); iter++) {

					vc_tmp = grid.gr.get_v_coord(*iter);
					tmp_pos = grid.get_pos_type(vc_tmp);

					if ((tmp_pos == B_EDGE) || (tmp_pos == L_B_CORNER)
							|| (tmp_pos == R_B_CORNER)){
						found = true;
						break;
					}
				} // for statement
				delete (test_mst_blue);
				if (found) {
					cout << "\n\nBLUE won the GAME  !!!\n\n";
					return hex_state::BLUE;
				}
			}
		} // outer for loop
	} // if statement ! full grid
	if (!full_grid)
		return hex_state::EMPTY;
	else
		return hex_state::BLUE;
}

//This is used in finding the max
template <class T>
bool vcMapTypeCompare (const T& x, const T& y){
	return (x.second < y.second);
}

/*
int evaluate_pos(Grid &gridx, vertex_coord& vc, hex_state player_color){
	list<vertex_coord> l_vc;
	int result = 0;
	if (player_color == hex_state::BLUE){
		// this goes up and down
		gridx.get_adjacent(vc, l_vc);

		for (auto iter = l_vc.begin(); iter != l_vc.end(); iter++) {
			int x = (*iter).first;
			int y = (*iter).second;
			//cout << g[x][y] << "  " << c << endl;
			if (gridx.g[x][y] == player_color) {
				if (abs(vc.second - y) == 1)
					// we are lengthing the chain
					result ++;
			}
		}
	}else {
		// RED goes from left to right
		gridx.get_adjacent(vc, l_vc);

		for (auto iter = l_vc.begin(); iter != l_vc.end(); iter++) {
			int x = (*iter).first;
			int y = (*iter).second;
			//cout << g[x][y] << "  " << c << endl;
			if (gridx.g[x][y] == player_color) {
				if (abs(vc.first - x) == 1)
					// we are lengthing the chain
					result ++;
			}
		}
	}
	return 0;
}

*/
vertex_coord monte_carlo(Grid & gridx, int simulations, hex_state player_color){
	vcMapType monte_result;
	vertex_coord start_vertex,vc_tmp,who_cares;
	//bool at_least_once = false;
	int computer_won = 0, opponent_won = 0;

	hex_state opp_color;
	if (player_color == hex_state::BLUE)
			opp_color = hex_state::RED;
	else
		opp_color = hex_state::BLUE;

	for (int i=0; i < gridx.get_dimension(); i++){
		for(int j = 0; j< gridx.get_dimension(); j++){
			start_vertex.first = i;
			start_vertex.second = j;
			if ( gridx.get_hex_state(start_vertex)== hex_state::EMPTY){


				for (int i=0; i < simulations; i++){
					Grid gridm(gridx);
					gridm.set(start_vertex,player_color);

					// fill_grid returns the first random position chosen for this player
					who_cares = fill_grid(gridm, opp_color);

					if (who_won(gridm) == player_color){
						monte_result[start_vertex]++;
						//at_least_once = true;
						computer_won++;
						//gridm.print_grid();
					}
					else
						opponent_won++;
				}
			}
		}
	}

	// debug print of the wins/trials
	for (auto i = monte_result.begin(); i != monte_result.end(); i++)
		cout << (*i).first.first << " " <<  (*i).first.second << " " << (*i).second << endl;

	// find the maximum
	auto pr = std::max_element(monte_result.begin(), monte_result.end(),vcMapTypeCompare<pair<vertex_coord,int> >);

	cout << "Max = " << pr->first.first << pr->first.second << endl;
	cout << "computer won " << computer_won << " opponent won " << opponent_won << endl;
	// return that location
	return pr->first;
}


int main() {
	srand((unsigned int) time(0));

	int size = 0;
	vector<vertex_coord> red_left;
	vector<vertex_coord> red_right;
	vector<vertex_coord> blue_top;
	vector<vertex_coord> blue_bottom;

	set<int> mst_red_set;
	set<int> mst_blue_set;

	//int comp_x, comp_y;
	int simulations = 1000;
	bool game_over = false;
	int num_v = 0;
	vertex_coord vc;
	 int x, y;
	char c = 'A';
	string mystr;
	cout << " \n\nLets Play a Game of Hex\n\n";
	cout << " \nEnter number of simulations for Monte-Carlo >";
	getline(cin,mystr);
	stringstream(mystr) >> simulations;
	//cin >> simulations;
	if ((simulations < 0) || (simulations > 1000)){
		simulations = 1000;
		cout << " Simulations set to 1000\n" << endl;
	}
	cout << " Please enter the size of the board [3,5,7,9,11,14 or 19] > ";


	//cin >> size;
	getline(cin,mystr);
	stringstream(mystr) >> size;

	if ((size < 3) || (size > 19)) {
		cout << " Invalid size for board- enter 7,9,11,14,or 19" << endl;
		return -1;
	}

	srand((unsigned int) time(0));

	// Initialize the grid
	Grid grid1(size);

	// Print Initial Grid
	grid1.print_grid();

	while (!game_over) {

		mst_red_set.clear();
		mst_blue_set.clear();
		if (num_v == size * size) {
			game_over = true;
			continue;
		}


		cout << " \n\nEnter you coordinates where (Row Col) eg. 1A or 3B > ";
		cin >> x >> c;
		y = static_cast<int>(toupper(c)) - 64;
		cout << y << endl;
		if ((x > size) || (x < 1) || (y < 1) || (y > size)) {
			cout << " Invalid Input " << "x = " << x << " y = " << y << endl;
			continue;
		}
		cout << " you chose (" << x << "," << static_cast<char>(toupper(c))
												<< ") or coord ("<<y-1<<","<<x-1<<")\n" << endl;

		// graph and grid indices start with 0 - so subtract 1
		vc.first = y - 1;
		vc.second = x - 1;

		// The following call - sets the grid location, locates adjacent nodes, appropriately adds the edges to
		// the graph.
		if (!grid1.set(vc, hex_state::RED)) {
			cout << "Grid Position (" << x << " , "
					<< static_cast<char>(toupper(c)) << ") already occupied \n";
			continue;
		}

		num_v++;

		//hex_pos_t red_pos = grid1.get_pos_type(vc);

		if (who_won(grid1,false) == hex_state::RED){
			grid1.print_grid();
			cout << "\n\nRED won the GAME - Computer (BLUE) Lost !!!\n\n";
			cout << " Total number of moves = " << num_v << endl;
			return 0;
		}


		if (num_v == size * size) {
			game_over = true;
			continue;
		}

		// This function finds the best next move for player of color
		vc = monte_carlo(grid1, simulations, hex_state::BLUE);

		if ((vc.first == -1) && (vc.second == -1)){
			grid1.print_grid();
			cout << "\n\nComputer (BLUE) Forfiets - No Path to win !!!\n\n";
			cout << " Total number of moves = " << num_v << endl;
			return 0;
		}
		cout << "BLUE computer chose " << vc.first << " "
				<< vc.second << endl;

		grid1.set(vc, hex_state::BLUE);


		num_v++;

		if (who_won(grid1,false) == hex_state::BLUE){
			grid1.print_grid();
			cout << "\n\nComputer (BLUE) won the GAME - You  Lost !!!\n\n";
			cout << " Total number of moves = " << num_v << endl;
			return 0;
		}

		grid1.print_grid();

	}
	grid1.print_grid();
	cout << "\nGAME OVER - No WINNERS - This is a bug, total moves " << num_v
			<< endl;
	;

}


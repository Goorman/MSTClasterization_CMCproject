#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <ctime>
#include <stdio.h>
#include <chrono>

using std::vector;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::string;

struct Node{
    double x;
    double y;
    int index;

    Node();
    Node(double _x, double _y, int _i) :
        x(_x), y(_y), index(_i){};
};

struct Edge{
    long v1;
    long v2;
    double distance;
    Edge(long _v1=-1, long _v2=-1, double _distance=0) : 
        v1(_v1), v2(_v2), distance(_distance){};
};

long find_set(long v, long* parent) {
    if (v == parent[v])
        return v;
    parent[v] = find_set(parent[v], parent);
    return parent[v];
}

void swap(long &a, long& b){
    long t = a;
    a = b;
    b = t;
}

void write_data(vector<Edge> edges, long nodes_number, const vector<Node>& nodes,
                string filename) {
    long parent[nodes_number];
    long rank[nodes_number];
    ofstream f(filename.c_str());
    for (long i = 0; i < nodes_number; ++i) {
        parent[i] = i;
        rank[i] = 0;
    }
    for (long i = 0; i < edges.size(); ++i) {
        long a = find_set(edges[i].v1, parent);
        long b = find_set(edges[i].v2, parent);
        if (a != b) {
            if (rank[a] < rank[b])
                swap(a, b);
            parent[b] = a;
            if (rank[a] == rank[b])
                ++rank[a];
        }
    }
    for (long i = 0; i < nodes_number; i++){
        find_set(i, parent);
    }
    for (int i = 0; i < nodes.size(); i++){
        f << nodes[i].x << "\t" << nodes[i].y << "\t" << parent[i] << endl;
    }
    f.close();
}

double get_distance(const Node& n1, const Node& n2){
    return sqrt(pow(n1.x - n2.x, 2) + pow(n1.y - n2.y, 2));
}

struct distcmp {
    bool operator() (const Edge& a, const Edge& b) {
        return a.distance < b.distance;
    }
};

void sort_tree(vector<Edge>& tree){
    std::sort(tree.begin(), tree.end(), distcmp());
}

Edge findMinimalOutgoingEdge(long comp_root, vector<Edge> edges,
                             long* comp_roots, long int N,
                             long comp_num) {
    Edge min_edge;
    if (comp_root != comp_roots[0]){
        min_edge = edges[comp_root * N + comp_roots[0]];
    }
    else{
        min_edge = edges[comp_root * N + comp_roots[1]];
    }
    for (long k = 0; k < comp_num; ++k) {
        if (comp_roots[k] == comp_root){
            continue;
        }
        if (min_edge.distance > edges[comp_root * N + comp_roots[k]].distance) {
            min_edge = edges[comp_root * N + comp_roots[k]];
        }
    }
    return min_edge;
}

long findComp(long* cc, long v1) {
    return cc[v1];
}

void connectComponentsWithEdge(long* cc, Edge edge, long N, vector<Edge>& edges) {
    long c1 = cc[edge.v1];
    long c2 = cc[edge.v2];
    #pragma omp parallel for
    for (long k = 0; k < N; ++k) {
        if (cc[k] == c2) {
            cc[k] = c1;
        }
        if (c2 == k){
            continue;
        }
        if (edges[c1 * N + k].distance > edges[c2 * N + k].distance) {
            edges[c1 * N + k] = edges[c2 * N + k];
            edges[k * N + c1] = edges[k * N + c2];
        }
    }
}

void get_clusters(vector<Node> nodes, string out_filename){
    cout << "size " << nodes.size() << endl;
    int clusters_number = 2;

    vector<Edge> tree;
    long int N = nodes.size();
    long cc[N];
    long comp_roots[N];
    long comp_num = N;
    vector<Edge> edges(N * N);
    #pragma omp parallel for
    for (long k = 0; k < N * N; ++k) {
        long i = k / N;
        long j = k % N;
        double distance = get_distance(nodes[i], nodes[j]);
        edges[i*N + j] = Edge(i, j, distance);
    }
    #pragma omp parallel for
    for (long i = 0; i < N; ++i) {
        cc[i] = i;
        comp_roots[i] = i;
    }

    while(comp_num > 1) {
        vector<Edge> mstedges(comp_num);
        #pragma omp parallel for
        for (int i = 0; i < comp_num; i++) {
            Edge e = findMinimalOutgoingEdge(comp_roots[i], edges, comp_roots,
                                             N, comp_num);
            mstedges[i] = e;
        }

        for (int i = 0; i < mstedges.size(); i++){
            if (findComp(cc, mstedges[i].v1) != findComp(cc, mstedges[i].v2)) {
                connectComponentsWithEdge(cc, mstedges[i], N, edges);
                tree.push_back(mstedges[i]);
            }
        }
        comp_num = 0;
        for (long i = 0; i < N; ++i) {
            if (cc[i] == i) {
            	comp_roots[comp_num] = i;
                ++comp_num;
            }
        }
        /*for (int i = 0; i < N; i++){
            cout << "cc [" <<i << "]=" << cc[i] << endl;
        }
        cout << "comp" << comp_num << endl;
        cout << "here2" << endl;*/
    }
    cout << "here1" << endl;
    sort_tree(tree);
    for (int i = 0; i < clusters_number-1; i++){
        tree.pop_back();
    }
    write_data(tree, N, nodes, out_filename);
}

void solvetask(string filename, string out_filename, int threadnum) {
    omp_set_num_threads(threadnum);

    ifstream f(filename.c_str());
    long n;
    f >> n;
    vector<Node> nodes;
    for (int i = 0; i < n; i++){
        double x, y;
        f >> x >> y;
        nodes.push_back(Node(x, y, i));
    }
    f.close();
    get_clusters(nodes, out_filename);
}

int main(int argc, const char * argv[]) {
    if (argc != 3){
        cout << "specify the in file name and out file name" << endl;
        return 1;
    }
    
    string filenamelist = argv[1];
    string out_filename = argv[2];
    int repeat_num = 5;

    int maxthreadnum = omp_get_max_threads(); 
    
    ifstream f(filenamelist.c_str()); 
    vector<string> filenames;
    int n = 0;
    f >> n;
    for (int i = 0; i < n; ++i) {
        string filename;
        f >> filename;
        filenames.push_back(filename);
    }
    f.close();
    ofstream out(out_filename.c_str());
    
    for (int a = 0; a < filenames.size(); ++a) {
        for (int b = 1; b <= maxthreadnum; ++b) {
	    std::chrono::time_point<std::chrono::system_clock> start, end;

	    start = std::chrono::system_clock::now();
	
            for (int i = 0; i < repeat_num; ++i){
            	solvetask(filenames[a], "tmp.txt", b);
            }
	    end = std::chrono::system_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - start;

	    out << a << " " << b << " " << elapsed_seconds.count()/repeat_num << std::endl;  
        }
    }
    out.close();

    return 0;
}

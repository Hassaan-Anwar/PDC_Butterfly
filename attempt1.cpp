// A minimal C++ adaptation of PARBUTTERFLY core algorithms
// Focus: Preprocessing, Wedge Retrieval, Vertex Butterfly Counting
// Updated to support Ligra-style general AdjacencyGraph format
// Modified to use C++11 compatible syntax

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <tuple>
#include <fstream>
#include <sstream>
#include <functional>
#include <map>
#include <chrono>  // For timing

using namespace std;
using namespace std::chrono;  // For timing
typedef int Vertex;
typedef pair<Vertex, Vertex> Edge;
typedef tuple<pair<Vertex, Vertex>, int, Vertex> Wedge;

// Hash function for pair<Vertex, Vertex>
struct pair_hash {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        return hash<T1>()(p.first) ^ (hash<T2>()(p.second) << 1);
    }
};

struct Graph {
    unordered_map<Vertex, vector<Vertex>> adj;
    vector<Vertex> vertices;
    unordered_map<Vertex, int> rank;
};

// Parse Ligra-style general adjacency graph format (AdjacencyGraph)
Graph loadLigraAdjacencyGraph(const string& filename) {
    auto start = high_resolution_clock::now();  // Declare start here
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        exit(1);
    }

    Graph G;
    string line;
    getline(file, line); // Skip the header "WeightedAdjacencyGraph"

    Vertex u, v;
    while (file >> u >> v) {
        G.adj[u].push_back(v);
    }

    unordered_set<Vertex> vertex_set;
    for (const auto& adj_pair : G.adj) {
        vertex_set.insert(adj_pair.first);
        for (Vertex neighbor : adj_pair.second) {
            vertex_set.insert(neighbor);
        }
    }
    G.vertices.assign(vertex_set.begin(), vertex_set.end());
    auto stop = high_resolution_clock::now();  // Stop timer
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Graph loaded in " << duration.count() << " milliseconds" << endl;

    return G;
}

// Ranking: sort all vertices based on a custom function (e.g., degree descending)
Graph preprocess(const Graph& G, function<bool(Vertex, Vertex)> rank_func) {
    Graph G_prime;
    G_prime.vertices = G.vertices;
    sort(G_prime.vertices.begin(), G_prime.vertices.end(), rank_func);
    for (size_t i = 0; i < G_prime.vertices.size(); ++i) {
        G_prime.rank[G_prime.vertices[i]] = i;
    }
    for (const auto& adj_pair : G.adj) {
        Vertex u = adj_pair.first;
        const vector<Vertex>& nbrs = adj_pair.second;
        vector<Vertex> sorted_nbrs = nbrs;
        sort(sorted_nbrs.begin(), sorted_nbrs.end(), [&](Vertex a, Vertex b) {
            return G_prime.rank[a] > G_prime.rank[b];
        });
        G_prime.adj[G_prime.rank[u]] = sorted_nbrs;
    }
    return G_prime;
}

vector<Wedge> getWedges(const Graph& G) {
    vector<Wedge> wedges;
    for (const auto& adj_pair : G.adj) {
        Vertex u1 = adj_pair.first;
        const vector<Vertex>& neighbors = adj_pair.second;
        for (Vertex v : neighbors) {
            if (G.adj.find(v) == G.adj.end()) continue;
            const vector<Vertex>& v_neighbors = G.adj.at(v);
            for (Vertex u2 : v_neighbors) {
                if (u1 != u2) {
                    wedges.emplace_back(make_pair(u1, u2), 1, v);
                }
            }
        }
    }
    return wedges;
}

int countTotalButterflies(const vector<Wedge>& wedges) {
    unordered_map<pair<Vertex, Vertex>, int, pair_hash> wedge_count;
    int total_butterflies = 0;

    for (const auto& wedge : wedges) {
        const pair<Vertex, Vertex>& endpoints = get<0>(wedge);
        wedge_count[endpoints]++;
    }

    for (const auto& wc_pair : wedge_count) {
        int d = wc_pair.second;
        total_butterflies += d * (d - 1);
    }

    unordered_map<Vertex, int> center_count;
    for (const auto& wedge : wedges) {
        const pair<Vertex, Vertex>& endpoints = get<0>(wedge);
        center_count[get<2>(wedge)] += wedge_count[endpoints] - 1;
    }
    for (const auto& cc_pair : center_count) {
        total_butterflies += cc_pair.second;
    }

    return total_butterflies;
}


int main() {
    
    auto total_start = high_resolution_clock::now();  // Declare total_start here

    auto load_start = high_resolution_clock::now();
    Graph G = loadLigraAdjacencyGraph("inputs/rMatGraph_J_5_500.txt");
    auto load_stop = high_resolution_clock::now();
    auto load_duration = duration_cast<milliseconds>(load_stop - load_start);
    cout << "Graph loading:      " << load_duration.count() << " ms\n";


    auto preprocess_start = high_resolution_clock::now();
    auto G_prime = preprocess(G, [&](Vertex a, Vertex b) {
        return G.adj[a].size() > G.adj[b].size();
    });
    auto preprocess_stop = high_resolution_clock::now();
    auto preprocess_duration = duration_cast<milliseconds>(preprocess_stop - preprocess_start);
    cout << "Preprocessing:      " << preprocess_duration.count() << " ms\n";


    auto wedges_start = high_resolution_clock::now();
    auto wedges = getWedges(G_prime);
    auto wedges_stop = high_resolution_clock::now();
    auto wedges_duration = duration_cast<milliseconds>(wedges_stop - wedges_start);
    cout << "Wedge computation:  " << wedges_duration.count() << " ms\n";


    auto count_start = high_resolution_clock::now();
    int total_butterflies = countTotalButterflies(wedges);
    auto count_stop = high_resolution_clock::now();
    auto count_duration = duration_cast<milliseconds>(count_stop - count_start);
    cout << "Butterfly counting: " << count_duration.count() << " ms\n";


    // Print total butterfly count instead of per vertex counts ...
    cout << "Total butterfly count: " << total_butterflies << endl;

    auto total_stop = high_resolution_clock::now();
    auto total_duration = duration_cast<milliseconds>(total_stop - total_start);

 /*   cout << "\nVertex butterfly counts:\n";
    for (const auto& bc_pair : butterfly_counts) {
        Vertex v = bc_pair.first;
        int count = bc_pair.second;
        cout << "Vertex " << v << ": " << count << endl;
    }
*/
    cout << "\nTiming Information:\n";
    cout << "--------------------------------\n";
    cout << "Graph loading:      " << load_duration.count() << " ms\n";
    cout << "Preprocessing:      " << preprocess_duration.count() << " ms\n";
    cout << "Wedge computation:  " << wedges_duration.count() << " ms\n";
    cout << "Butterfly counting: " << count_duration.count() << " ms\n";
    cout << "--------------------------------\n";
    cout << "Total time:         " << total_duration.count() << " ms\n";

    return 0;
}
// A minimal C++ adaptation of PARBUTTERFLY core algorithms
// Focus: Preprocessing, Wedge Retrieval, Vertex Butterfly Counting
// Updated to support Ligra-style general AdjacencyGraph format
// Modified to use C++11 compatible syntax and OpenMP for parallelism

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
#include <chrono>
#include <omp.h>  // OpenMP

using namespace std;
using namespace std::chrono;
typedef int Vertex;
typedef pair<Vertex, Vertex> Edge;
typedef tuple<pair<Vertex, Vertex>, int, Vertex> Wedge;

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

Graph loadLigraAdjacencyGraph(const string& filename) {
    auto start = high_resolution_clock::now();
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        exit(1);
    }

    Graph G;
    string line;
    getline(file, line); // Skip the header

    Vertex u, v;
    while (file >> u >> v) {
        #pragma omp critical
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

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Graph loaded in " << duration.count() << " milliseconds" << endl;

    return G;
}

Graph preprocess(const Graph& G, function<bool(Vertex, Vertex)> rank_func) {
    Graph G_prime;
    G_prime.vertices = G.vertices;
    sort(G_prime.vertices.begin(), G_prime.vertices.end(), rank_func);
    for (size_t i = 0; i < G_prime.vertices.size(); ++i) {
        G_prime.rank[G_prime.vertices[i]] = i;
    }
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < G.vertices.size(); ++i) {
        Vertex u = G.vertices[i];
        if (!G.adj.count(u)) continue;
        const auto& nbrs = G.adj.at(u);
        vector<Vertex> sorted_nbrs = nbrs;
        sort(sorted_nbrs.begin(), sorted_nbrs.end(), [&](Vertex a, Vertex b) {
            int rank_a = G_prime.rank.count(a) ? G_prime.rank[a] : -1;
            int rank_b = G_prime.rank.count(b) ? G_prime.rank[b] : -1;
            return rank_a > rank_b;
        });
        #pragma omp critical
        G_prime.adj[G_prime.rank[u]] = sorted_nbrs;
    }
    return G_prime;
}

vector<Wedge> getWedges(const Graph& G) {
    vector<Wedge> wedges;
    #pragma omp parallel
    {
        vector<Wedge> local_wedges;
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < G.vertices.size(); ++i) {
            Vertex u1 = G.vertices[i];
            if (!G.adj.count(u1)) continue;
            const auto& neighbors = G.adj.at(u1);
            for (Vertex v : neighbors) {
                if (!G.adj.count(v)) continue;
                const auto& v_neighbors = G.adj.at(v);
                for (Vertex u2 : v_neighbors) {
                    if (u1 != u2) {
                        local_wedges.emplace_back(make_pair(u1, u2), 1, v);
                    }
                }
            }
        }
        #pragma omp critical
        wedges.insert(wedges.end(), local_wedges.begin(), local_wedges.end());
    }
    return wedges;
}

int countTotalButterflies(const vector<Wedge>& wedges) {
    unordered_map<pair<Vertex, Vertex>, int, pair_hash> wedge_count;
    int total_butterflies = 0;

    #pragma omp parallel
    {
        unordered_map<pair<Vertex, Vertex>, int, pair_hash> local_count;
        #pragma omp for nowait
        for (int i = 0; i < wedges.size(); ++i) {
            const auto& endpoints = get<0>(wedges[i]);
            local_count[endpoints]++;
        }
        #pragma omp critical
        for (auto it = local_count.begin(); it != local_count.end(); ++it) {
            wedge_count[it->first] += it->second;
        }
    }

    #pragma omp parallel for reduction(+:total_butterflies)
    for (int i = 0; i < wedges.size(); ++i) {
        const auto& endpoints = get<0>(wedges[i]);
        int d = wedge_count[endpoints];
        total_butterflies += d * (d - 1);
    }

    unordered_map<Vertex, int> center_count;
    #pragma omp parallel
    {
        unordered_map<Vertex, int> local_center;
        #pragma omp for nowait
        for (int i = 0; i < wedges.size(); ++i) {
            const auto& endpoints = get<0>(wedges[i]);
            Vertex center = get<2>(wedges[i]);
            local_center[center] += wedge_count[endpoints] - 1;
        }
        #pragma omp critical
        for (auto it = local_center.begin(); it != local_center.end(); ++it) {
            center_count[it->first] += it->second;
        }
    }

    for (auto it = center_count.begin(); it != center_count.end(); ++it) {
        total_butterflies += it->second;
    }

    return total_butterflies;
}

int main() {
    auto total_start = high_resolution_clock::now();

    auto load_start = high_resolution_clock::now();
    Graph G = loadLigraAdjacencyGraph("inputs/rMatGraph_J_5_500.txt");
    auto load_stop = high_resolution_clock::now();
    auto load_duration = duration_cast<milliseconds>(load_stop - load_start);

    auto preprocess_start = high_resolution_clock::now();
    auto G_prime = preprocess(G, [&](Vertex a, Vertex b) {
        size_t deg_a = G.adj.count(a) ? G.adj.at(a).size() : 0;
        size_t deg_b = G.adj.count(b) ? G.adj.at(b).size() : 0;
        return deg_a > deg_b;
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

    cout << "Total butterfly count: " << total_butterflies << endl;

    auto total_stop = high_resolution_clock::now();
    auto total_duration = duration_cast<milliseconds>(total_stop - total_start);

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

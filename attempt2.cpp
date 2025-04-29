// A minimal C++ adaptation of PARBUTTERFLY core algorithms
// Focus: Preprocessing, Wedge Retrieval, Vertex Butterfly Counting
// Updated to support Ligra-style bipartite input format

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

using namespace std;

typedef int Vertex;
typedef pair<Vertex, Vertex> Edge;
typedef tuple<pair<Vertex, Vertex>, int, Vertex> Wedge;

struct Graph {
    unordered_map<Vertex, vector<Vertex>> adj;
    vector<Vertex> vertices;
    unordered_map<Vertex, int> rank;
};

// Parse Ligra-style bipartite format
Graph loadLigraBipartiteGraph(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        exit(1);
    }

    int nv, mv, nu, mu;
    file >> nv >> mv >> nu >> mu;

    vector<int> offsetV(nv), offsetU(nu);
    vector<int> edgesV(mv), edgesU(mu);

    for (int i = 0; i < nv; ++i) file >> offsetV[i];
    for (int i = 0; i < mv; ++i) file >> edgesV[i];
    for (int i = 0; i < nu; ++i) file >> offsetU[i];
    for (int i = 0; i < mu; ++i) file >> edgesU[i];

    Graph G;
    // Build adj from V to U
    for (int i = 0; i < nv; ++i) {
        int start = offsetV[i];
        int end = (i + 1 < nv) ? offsetV[i + 1] : mv;
        for (int j = start; j < end; ++j) {
            G.adj[i].push_back(edgesV[j]);
        }
        G.vertices.push_back(i);
    }

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
    for (auto& [u, nbrs] : G.adj) {
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
    for (const auto& [u1, neighbors] : G.adj) {
        for (Vertex v : neighbors) {
            if (G.adj.find(v) == G.adj.end()) continue;
            for (Vertex u2 : G.adj.at(v)) {
                if (u1 != u2) {
                    wedges.emplace_back(make_pair(u1, u2), 1, v);
                }
            }
        }
    }
    return wedges;
}

unordered_map<Vertex, int> countVertexButterflies(const vector<Wedge>& wedges) {
    unordered_map<pair<Vertex, Vertex>, int, hash<pair<int,int>>> wedge_count;
    unordered_map<Vertex, int> butterfly_count;

    for (const auto& [endpoints, freq, center] : wedges) {
        wedge_count[endpoints]++;
    }

    for (const auto& [endpoints, d] : wedge_count) {
        Vertex u1 = endpoints.first;
        Vertex u2 = endpoints.second;
        butterfly_count[u1] += d * (d - 1) / 2;
        butterfly_count[u2] += d * (d - 1) / 2;
    }

    unordered_map<Vertex, int> center_count;
    for (const auto& [endpoints, freq, center] : wedges) {
        center_count[center] += wedge_count[endpoints] - 1;
    }
    for (const auto& [v, c] : center_count) {
        butterfly_count[v] += c;
    }

    return butterfly_count;
}

int main() {
    Graph G = loadLigraBipartiteGraph("inputs/rMatGraph_J_5_100.txt");

    // Degree-based ranking
    auto G_prime = preprocess(G, [&](Vertex a, Vertex b) {
        return G.adj[a].size() > G.adj[b].size();
    });

    auto wedges = getWedges(G_prime);
    auto butterfly_counts = countVertexButterflies(wedges);

    cout << "Vertex butterfly counts:\n";
    for (auto& [v, count] : butterfly_counts) {
        cout << "Vertex " << v << ": " << count << endl;
    }
    return 0;
}

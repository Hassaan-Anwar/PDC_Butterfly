#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <tuple>

using namespace std;

typedef pair<int, int> Edge;
typedef unordered_map<int, unordered_set<int>> BipartiteGraph;

// Load bipartite edges from file (u v)
BipartiteGraph loadGraph(const string& filename, double sampleProb) {
    ifstream infile(filename);
    string line;
    BipartiteGraph graph;
    mt19937 rng(random_device{}());
    uniform_real_distribution<double> dist(0.0, 1.0);

    while (getline(infile, line)) {
        istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) continue;

        // Sparsify: sample edge
        if (dist(rng) <= sampleProb) {
            graph[u].insert(v);  // u ∈ U, v ∈ V
        }
    }
    return graph;
}

// Count butterflies in the (sparse) bipartite graph
int countButterflies(const BipartiteGraph& graph) {
    unordered_map<int, vector<int>> common;

    for (const auto& [u1, neighbors1] : graph) {
        for (const auto& [u2, neighbors2] : graph) {
            if (u1 >= u2) continue;

            // Find common neighbors between u1 and u2 (in V)
            vector<int> intersection;
            for (int v : neighbors1) {
                if (neighbors2.find(v) != neighbors2.end()) {
                    intersection.push_back(v);
                }
            }

            int k = intersection.size();
            if (k >= 2) {
                // Number of butterflies is C(k, 2)
                common[u1].push_back(k * (k - 1) / 2);
            }
        }
    }

    // Sum up
    int total = 0;
    for (auto& [u, counts] : common) {
        for (int c : counts)
            total += c;
    }

    return total;
}

int main() {
    string filename = "graph.txt";  // Your bipartite edge list file
    double p = 0.1;  // Sampling probability

    BipartiteGraph sparseGraph = loadGraph(filename, p);
    int sparseCount = countButterflies(sparseGraph);

    // Scale up the count (approximate)
    double estimatedCount = sparseCount / pow(p, 4);  // 4 edges per butterfly

    cout << "Approximate number of butterflies: " << round(estimatedCount) << endl;

    return 0;
}

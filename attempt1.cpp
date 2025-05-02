<<<<<<< HEAD
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
=======
// Optimized C++ implementation of PARBUTTERFLY using OpenMPI
// Focus: Preprocessing, Wedge Retrieval, Vertex Butterfly Counting
// Parallelized for distributed computing environments

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
#include <mpi.h>   // OpenMPI

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
Graph loadLigraAdjacencyGraph(const string& filename, int rank, int size) {
    auto start = high_resolution_clock::now();
    
    // Only rank 0 reads the file and broadcasts the graph
    Graph G;
    
    if (rank == 0) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file!" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

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
    }
    
    // Serialize and broadcast the graph structure
    // First broadcast number of vertices and edges for allocation
    int num_vertices = (rank == 0) ? G.vertices.size() : 0;
    int num_adj_entries = 0;
    
    if (rank == 0) {
        for (const auto& adj_pair : G.adj) {
            num_adj_entries += adj_pair.second.size() + 1; // +1 for the vertex itself
        }
    }
    
    MPI_Bcast(&num_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_adj_entries, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        // Allocate space for vertices
        G.vertices.resize(num_vertices);
    }
    
    // Broadcast vertices
    MPI_Bcast(G.vertices.data(), num_vertices, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Serialize and broadcast adjacency lists
    vector<int> serialized_graph;
    if (rank == 0) {
        serialized_graph.reserve(num_adj_entries);
        for (const auto& adj_pair : G.adj) {
            serialized_graph.push_back(adj_pair.first);
            serialized_graph.push_back(adj_pair.second.size());
            for (Vertex v : adj_pair.second) {
                serialized_graph.push_back(v);
            }
        }
    } else {
        serialized_graph.resize(num_adj_entries);
    }
    
    MPI_Bcast(serialized_graph.data(), num_adj_entries, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Deserialize if not rank 0
    if (rank != 0) {
        G.adj.clear();
        size_t idx = 0;
        while (idx < serialized_graph.size()) {
            Vertex u = serialized_graph[idx++];
            int degree = serialized_graph[idx++];
            vector<Vertex> neighbors;
            neighbors.reserve(degree);
            for (int i = 0; i < degree; i++) {
                neighbors.push_back(serialized_graph[idx++]);
            }
            G.adj[u] = neighbors;
        }
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    if (rank == 0) {
        cout << "Graph loaded and distributed in " << duration.count() << " milliseconds" << endl;
    }

    return G;
}

// Ranking: sort all vertices based on a custom function (e.g., degree descending)
Graph preprocess(const Graph& G, function<bool(Vertex, Vertex)> rank_func, int rank, int size) {
    Graph G_prime;
    G_prime.vertices = G.vertices;
    
    // Sort vertices - each process sorts its own copy
    sort(G_prime.vertices.begin(), G_prime.vertices.end(), rank_func);
    
    // Create rank map
    for (size_t i = 0; i < G_prime.vertices.size(); ++i) {
        G_prime.rank[G_prime.vertices[i]] = i;
    }
    
    // Chunk the adjacency lists for parallel processing
    int chunk_size = G.adj.size() / size;
    int start_idx = rank * chunk_size;
    int end_idx = (rank == size - 1) ? G.adj.size() : start_idx + chunk_size;
    
    // Process assigned adjacency lists
    int counter = 0;
    for (const auto& adj_pair : G.adj) {
        if (counter >= start_idx && counter < end_idx) {
            Vertex u = adj_pair.first;
            const vector<Vertex>& nbrs = adj_pair.second;
            vector<Vertex> sorted_nbrs = nbrs;
            sort(sorted_nbrs.begin(), sorted_nbrs.end(), [&](Vertex a, Vertex b) {
                return G_prime.rank[a] > G_prime.rank[b];
            });
            G_prime.adj[G_prime.rank[u]] = sorted_nbrs;
        }
        counter++;
    }
    
    // Gather all processed adjacency lists
    // This is a simplified approach - in a real implementation, you would use MPI structures
    // to efficiently gather the distributed adjacency lists
    
    // For each process, serialize its chunk of adjacency lists
    int local_adj_size = 0;
    for (const auto& adj_pair : G_prime.adj) {
        local_adj_size += 2 + adj_pair.second.size(); // vertex + size + neighbors
    }
    
    vector<int> local_serialized(local_adj_size);
    int idx = 0;
    for (const auto& adj_pair : G_prime.adj) {
        local_serialized[idx++] = adj_pair.first;
        local_serialized[idx++] = adj_pair.second.size();
        for (Vertex v : adj_pair.second) {
            local_serialized[idx++] = v;
        }
    }
    
    // Gather sizes from all processes
    vector<int> all_sizes(size);
    MPI_Allgather(&local_adj_size, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    // Calculate displacements
    vector<int> displacements(size);
    displacements[0] = 0;
    for (int i = 1; i < size; i++) {
        displacements[i] = displacements[i-1] + all_sizes[i-1];
    }
    
    // Calculate total size
    int total_size = 0;
    for (int s : all_sizes) {
        total_size += s;
    }
    
    // Gather all serialized data
    vector<int> global_serialized(total_size);
    MPI_Allgatherv(local_serialized.data(), local_adj_size, MPI_INT,
                  global_serialized.data(), all_sizes.data(), displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Clear current adjacency lists and rebuild from gathered data
    G_prime.adj.clear();
    idx = 0;
    while (idx < global_serialized.size()) {
        Vertex u = global_serialized[idx++];
        int degree = global_serialized[idx++];
        vector<Vertex> neighbors;
        neighbors.reserve(degree);
        for (int i = 0; i < degree; i++) {
            neighbors.push_back(global_serialized[idx++]);
        }
        G_prime.adj[u] = neighbors;
    }
    
    return G_prime;
}

vector<Wedge> getWedges(const Graph& G, int rank, int size) {
    vector<Wedge> local_wedges;
    
    // Assign vertices to processes in a round-robin fashion
    for (size_t i = rank; i < G.vertices.size(); i += size) {
        Vertex u1 = G.vertices[i];
        
        // Skip if u1 is not in the adjacency list
        if (G.adj.find(u1) == G.adj.end()) continue;
        
        const vector<Vertex>& neighbors = G.adj.at(u1);
        for (Vertex v : neighbors) {
            if (G.adj.find(v) == G.adj.end()) continue;
            const vector<Vertex>& v_neighbors = G.adj.at(v);
            for (Vertex u2 : v_neighbors) {
                if (u1 != u2) {
                    local_wedges.emplace_back(make_pair(u1, u2), 1, v);
                }
            }
        }
    }
    
    // Gather wedge counts from all processes
    int local_wedge_count = local_wedges.size();
    vector<int> all_wedge_counts(size);
    
    MPI_Allgather(&local_wedge_count, 1, MPI_INT, all_wedge_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    // Calculate total wedges and displacements
    vector<int> displacements(size);
    displacements[0] = 0;
    for (int i = 1; i < size; i++) {
        displacements[i] = displacements[i-1] + all_wedge_counts[i-1];
    }
    
    int total_wedges = 0;
    for (int count : all_wedge_counts) {
        total_wedges += count;
    }
    
    // Serialize local wedges for gathering
    // Format: u1, u2, dummy (always 1), v
    vector<int> serialized_local_wedges(local_wedge_count * 4);
    for (size_t i = 0; i < local_wedges.size(); i++) {
        const Wedge& w = local_wedges[i];
        const pair<Vertex, Vertex>& endpoints = get<0>(w);
        int dummy = get<1>(w);
        Vertex center = get<2>(w);
        
        serialized_local_wedges[i*4] = endpoints.first;
        serialized_local_wedges[i*4+1] = endpoints.second;
        serialized_local_wedges[i*4+2] = dummy;
        serialized_local_wedges[i*4+3] = center;
    }
    
    // Gather all wedges
    vector<int> serialized_all_wedges(total_wedges * 4);
    vector<int> recv_counts(size);
    for (int i = 0; i < size; i++) {
        recv_counts[i] = all_wedge_counts[i] * 4;
    }
    vector<int> recv_displacements(size);
    for (int i = 0; i < size; i++) {
        recv_displacements[i] = displacements[i] * 4;
    }
    
    MPI_Allgatherv(serialized_local_wedges.data(), local_wedge_count * 4, MPI_INT,
                  serialized_all_wedges.data(), recv_counts.data(), recv_displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Deserialize to vector<Wedge>
    vector<Wedge> all_wedges;
    all_wedges.reserve(total_wedges);
    
    for (int i = 0; i < total_wedges; i++) {
        Vertex u1 = serialized_all_wedges[i*4];
        Vertex u2 = serialized_all_wedges[i*4+1];
        int dummy = serialized_all_wedges[i*4+2];
        Vertex center = serialized_all_wedges[i*4+3];
        
        all_wedges.emplace_back(make_pair(u1, u2), dummy, center);
    }
    
    return all_wedges;
}

unordered_map<Vertex, int> countVertexButterflies(const vector<Wedge>& wedges, int rank, int size) {
    unordered_map<pair<Vertex, Vertex>, int, pair_hash> local_wedge_count;
    unordered_map<Vertex, int> local_butterfly_count;
    
    // Distribute wedges among processes
    int wedges_per_proc = wedges.size() / size;
    int start_idx = rank * wedges_per_proc;
    int end_idx = (rank == size - 1) ? wedges.size() : start_idx + wedges_per_proc;
    
    // Count wedges between vertex pairs
    for (int i = start_idx; i < end_idx; i++) {
        const auto& wedge = wedges[i];
        const pair<Vertex, Vertex>& endpoints = get<0>(wedge);
        local_wedge_count[endpoints]++;
    }
    
    // Gather all wedge counts
    // Serialize local_wedge_count
    vector<int> serialized_local_counts;
    serialized_local_counts.reserve(local_wedge_count.size() * 3); // u1, u2, count
    
    for (const auto& wc_pair : local_wedge_count) {
        const pair<Vertex, Vertex>& endpoints = wc_pair.first;
        int count = wc_pair.second;
        
        serialized_local_counts.push_back(endpoints.first);
        serialized_local_counts.push_back(endpoints.second);
        serialized_local_counts.push_back(count);
    }
    
    // Gather sizes of serialized data
    int local_size = serialized_local_counts.size();
    vector<int> all_sizes(size);
    
    MPI_Allgather(&local_size, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    // Calculate displacements
    vector<int> displacements(size);
    displacements[0] = 0;
    for (int i = 1; i < size; i++) {
        displacements[i] = displacements[i-1] + all_sizes[i-1];
    }
    
    // Calculate total size
    int total_size = 0;
    for (int s : all_sizes) {
        total_size += s;
    }
    
    // Gather all serialized data
    vector<int> all_serialized_counts(total_size);
    
    MPI_Allgatherv(serialized_local_counts.data(), local_size, MPI_INT,
                  all_serialized_counts.data(), all_sizes.data(), displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Deserialize and consolidate wedge counts
    unordered_map<pair<Vertex, Vertex>, int, pair_hash> wedge_count;
    
    for (size_t i = 0; i < all_serialized_counts.size(); i += 3) {
        Vertex u1 = all_serialized_counts[i];
        Vertex u2 = all_serialized_counts[i+1];
        int count = all_serialized_counts[i+2];
        
        wedge_count[make_pair(u1, u2)] += count;
    }
    
    // Calculate butterfly counts from wedge counts
    // Assign vertex pairs to processes
    int pairs_per_proc = wedge_count.size() / size;
    
    int pair_idx = 0;
    for (const auto& wc_pair : wedge_count) {
        if (pair_idx % size == rank) {
            const pair<Vertex, Vertex>& endpoints = wc_pair.first;
            int d = wc_pair.second;
            Vertex u1 = endpoints.first;
            Vertex u2 = endpoints.second;
            
            // Each endpoint gets d*(d-1)/2 butterflies
            local_butterfly_count[u1] += d * (d - 1) / 2;
            local_butterfly_count[u2] += d * (d - 1) / 2;
        }
        pair_idx++;
    }
    
    // Count butterflies for center vertices
    unordered_map<Vertex, int> local_center_count;
    
    // Divide wedges among processes
    for (int i = start_idx; i < end_idx; i++) {
        const auto& wedge = wedges[i];
        const pair<Vertex, Vertex>& endpoints = get<0>(wedge);
        Vertex center = get<2>(wedge);
        local_center_count[center] += wedge_count[endpoints] - 1;
    }
    
    // Add center butterfly counts
    for (const auto& cc_pair : local_center_count) {
        Vertex v = cc_pair.first;
        int c = cc_pair.second;
        local_butterfly_count[v] += c;
    }
    
    // Gather all butterfly counts
    // Serialize local butterfly counts
    vector<int> serialized_local_bf_counts;
    serialized_local_bf_counts.reserve(local_butterfly_count.size() * 2); // vertex, count
    
    for (const auto& bc_pair : local_butterfly_count) {
        Vertex v = bc_pair.first;
        int count = bc_pair.second;
        
        serialized_local_bf_counts.push_back(v);
        serialized_local_bf_counts.push_back(count);
    }
    
    // Gather sizes of serialized data
    local_size = serialized_local_bf_counts.size();
    all_sizes.resize(size);
    
    MPI_Allgather(&local_size, 1, MPI_INT, all_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    // Calculate displacements
    displacements.resize(size);
    displacements[0] = 0;
    for (int i = 1; i < size; i++) {
        displacements[i] = displacements[i-1] + all_sizes[i-1];
    }
    
    // Calculate total size
    total_size = 0;
    for (int s : all_sizes) {
        total_size += s;
    }
    
    // Gather all serialized data
    vector<int> all_serialized_bf_counts(total_size);
    
    MPI_Allgatherv(serialized_local_bf_counts.data(), local_size, MPI_INT,
                  all_serialized_bf_counts.data(), all_sizes.data(), displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Deserialize and consolidate butterfly counts
    unordered_map<Vertex, int> butterfly_count;
    
    for (size_t i = 0; i < all_serialized_bf_counts.size(); i += 2) {
        Vertex v = all_serialized_bf_counts[i];
        int count = all_serialized_bf_counts[i+1];
        
        butterfly_count[v] += count;
    }
    
    return butterfly_count;
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        cout << "Running with " << size << " MPI processes" << endl;
    }
    
    auto total_start = high_resolution_clock::now();
    
    // Set default input file or take from command line
    string input_file = "inputs/rMatGraph_J_5_500.txt";
    if (argc > 1) {
        input_file = argv[1];
    }
    
    auto load_start = high_resolution_clock::now();
    Graph G = loadLigraAdjacencyGraph(input_file, rank, size);
    auto load_stop = high_resolution_clock::now();
    auto load_duration = duration_cast<milliseconds>(load_stop - load_start);
    
    auto preprocess_start = high_resolution_clock::now();
    auto G_prime = preprocess(G, [&](Vertex a, Vertex b) {
        size_t size_a = G.adj.count(a) ? G.adj[a].size() : 0;
        size_t size_b = G.adj.count(b) ? G.adj[b].size() : 0;
        return size_a > size_b;
    }, rank, size);
    auto preprocess_stop = high_resolution_clock::now();
    auto preprocess_duration = duration_cast<milliseconds>(preprocess_stop - preprocess_start);
    
    auto wedges_start = high_resolution_clock::now();
    auto wedges = getWedges(G_prime, rank, size);
    auto wedges_stop = high_resolution_clock::now();
    auto wedges_duration = duration_cast<milliseconds>(wedges_stop - wedges_start);
    
    auto count_start = high_resolution_clock::now();
    auto butterfly_counts = countVertexButterflies(wedges, rank, size);
    auto count_stop = high_resolution_clock::now();
    auto count_duration = duration_cast<milliseconds>(count_stop - count_start);
    
    auto total_stop = high_resolution_clock::now();
    auto total_duration = duration_cast<milliseconds>(total_stop - total_start);
    
    // Calculate total butterfly count
    long long local_total = 0;
    for (const auto& bc_pair : butterfly_counts) {
        local_total += bc_pair.second;
    }
    
    long long global_total = 0;
    MPI_Reduce(&local_total, &global_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    // Since each butterfly is counted 4 times (once for each vertex), divide by 4
    global_total /= 4;
    
    if (rank == 0) {
        cout << "\nVertex butterfly counts:\n";
        for (const auto& bc_pair : butterfly_counts) {
            Vertex v = bc_pair.first;
            int count = bc_pair.second;
            cout << "Vertex " << v << ": " << count << endl;
        }
        
        cout << "\nSummary:\n";
        cout << "--------------------------------\n";
        cout << "Total wedges:       " << wedges.size() << endl;
        cout << "Total butterflies:  " << global_total << endl;
        cout << "--------------------------------\n";
        
        cout << "\nTiming Information:\n";
        cout << "--------------------------------\n";
        cout << "Graph loading:      " << load_duration.count() << " ms\n";
        cout << "Preprocessing:      " << preprocess_duration.count() << " ms\n";
        cout << "Wedge computation:  " << wedges_duration.count() << " ms\n";
        cout << "Butterfly counting: " << count_duration.count() << " ms\n";
        cout << "--------------------------------\n";
        cout << "Total time:         " << total_duration.count() << " ms\n";
    }
    
    // Finalize MPI
    MPI_Finalize();
    
    return 0;
>>>>>>> e9fe162 (Refactor and enhance MPI-based butterfly counting algorithm)
}
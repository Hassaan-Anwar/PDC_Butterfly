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
#include <mpi.h>   // Include MPI header
#include <metis.h>
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
Graph loadLigraAdjacencyGraph(const string& filename, int rank, int world_size) {
    auto start = high_resolution_clock::now();
    
    Graph G;
    
    // Only process 0 reads the file and broadcasts the data
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
    
    // Serialize the graph for broadcasting
    vector<int> serialized_data;
    
    if (rank == 0) {
        // First element is the number of vertices
        serialized_data.push_back(G.vertices.size());
        
        // Add all vertices
        serialized_data.insert(serialized_data.end(), G.vertices.begin(), G.vertices.end());
        
        // For each vertex, add adjacency list
        for (const auto& v : G.vertices) {
            if (G.adj.find(v) != G.adj.end()) {
                // Add number of neighbors
                serialized_data.push_back(G.adj[v].size());
                // Add all neighbors
                serialized_data.insert(serialized_data.end(), G.adj[v].begin(), G.adj[v].end());
            } else {
                // No neighbors
                serialized_data.push_back(0);
            }
        }
    }
    
    // Broadcast the size of serialized data first
    int data_size = serialized_data.size();
    MPI_Bcast(&data_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Resize serialized_data on non-root processes
    if (rank != 0) {
        serialized_data.resize(data_size);
    }
    
    // Broadcast the serialized data
    MPI_Bcast(serialized_data.data(), data_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Deserialize the data on all processes
    if (rank != 0) {
        int pos = 0;
        int num_vertices = serialized_data[pos++];
        
        // Get all vertices
        G.vertices.resize(num_vertices);
        for (int i = 0; i < num_vertices; i++) {
            G.vertices[i] = serialized_data[pos++];
        }
        
        // Get adjacency lists
        for (int i = 0; i < num_vertices; i++) {
            Vertex v = G.vertices[i];
            int num_neighbors = serialized_data[pos++];
            
            // Add all neighbors
            for (int j = 0; j < num_neighbors; j++) {
                G.adj[v].push_back(serialized_data[pos++]);
            }
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
        G_prime.adj[u] = sorted_nbrs;
    }
    return G_prime;
}

vector<Wedge> getWedges(const Graph& G, int rank, int world_size, 
    const unordered_map<Vertex, int>& vertex_partition) {
vector<Wedge> wedges;

// Collect vertices assigned to this partition via METIS
vector<Vertex> my_vertices;
for (Vertex u : G.vertices) {
if (vertex_partition.at(u) == rank) {
my_vertices.push_back(u);
}
}

// Process assigned vertices using FULL adjacency lists
#pragma omp parallel
{
    vector<Wedge> local_wedges;
    #pragma omp for nowait
    for (int i = 0; i < my_vertices.size(); ++i) {
        Vertex u1 = my_vertices[i];
        if (!G.adj.count(u1)) continue;

        const auto& u1_neighbors = G.adj.at(u1);
        for (Vertex v : u1_neighbors) {
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

    // Gather wedge counts
    int local_wedge_count = wedges.size();
    vector<int> all_wedge_counts(world_size);
    MPI_Allgather(&local_wedge_count, 1, MPI_INT, all_wedge_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    vector<int> displacements(world_size, 0);
    for (int i = 1; i < world_size; i++) {
        displacements[i] = displacements[i-1] + all_wedge_counts[i-1];
    }
    
    // Create a serialized format for wedges
    vector<int> serialized_wedges;
    for (const auto& wedge : wedges) {
        const pair<Vertex, Vertex>& endpoints = get<0>(wedge);
        serialized_wedges.push_back(endpoints.first);
        serialized_wedges.push_back(endpoints.second);
        serialized_wedges.push_back(get<1>(wedge));
        serialized_wedges.push_back(get<2>(wedge));
    }
    
    // Calculate total wedges
    int total_wedges = 0;
    for (int count : all_wedge_counts) {
        total_wedges += count;
    }
    
    vector<int> all_serialized_wedges(total_wedges * 4);
    vector<int> recv_counts(world_size);
    for (int i = 0; i < world_size; i++) {
        recv_counts[i] = all_wedge_counts[i] * 4;
    }
    vector<int> send_displacements(world_size);
    for (int i = 1; i < world_size; i++) {
        send_displacements[i] = send_displacements[i-1] + recv_counts[i-1];
    }
    
    // All-gather the serialized wedges
    MPI_Allgatherv(serialized_wedges.data(), serialized_wedges.size(), MPI_INT,
                  all_serialized_wedges.data(), recv_counts.data(), send_displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Deserialize wedges
    vector<Wedge> all_wedges;
    for (int i = 0; i < all_serialized_wedges.size(); i += 4) {
        Vertex u1 = all_serialized_wedges[i];
        Vertex u2 = all_serialized_wedges[i+1];
        int weight = all_serialized_wedges[i+2];
        Vertex center = all_serialized_wedges[i+3];
        all_wedges.emplace_back(make_pair(u1, u2), weight, center);
    }
    
    return all_wedges;
}

unordered_map<Vertex, int> countVertexButterflies(const vector<Wedge>& wedges, int rank, int world_size, 
                                          const unordered_map<Vertex, int>& vertex_partition) {
    unordered_map<pair<Vertex, Vertex>, int, pair_hash> wedge_count;
    unordered_map<Vertex, int> butterfly_count;
    
    // Filter wedges according to the current partition
    vector<Wedge> my_wedges;
    #pragma omp parallel
    {
        vector<Wedge> local_my_wedges;
        #pragma omp for nowait
        for (int i = 0; i < wedges.size(); ++i) {
            const auto& wedge = wedges[i];
            const auto& endpoints = get<0>(wedge);
            Vertex u1 = endpoints.first;
            Vertex u2 = endpoints.second;
            Vertex center = get<2>(wedge);
            if (vertex_partition.at(u1) == rank ||
                vertex_partition.at(u2) == rank ||
                vertex_partition.at(center) == rank) {
                local_my_wedges.push_back(wedge);
            }
        }
        #pragma omp critical
        my_wedges.insert(my_wedges.end(), local_my_wedges.begin(), local_my_wedges.end());
    }
    
    // Local wedge counting using partition-filtered wedges
    #pragma omp parallel
    {
        unordered_map<pair<Vertex, Vertex>, int, pair_hash> local_wedge_count;
        #pragma omp for nowait
        for (int i = 0; i < my_wedges.size(); ++i) {
            const auto& wedge = my_wedges[i];
            const auto& endpoints = get<0>(wedge);
            local_wedge_count[endpoints]++;
        }
        #pragma omp critical
        for (const auto& p : local_wedge_count) {
            wedge_count[p.first] += p.second;
        }
    }
    
    
    // All-reduce wedge counts
    vector<pair<pair<Vertex, Vertex>, int>> local_wedge_count_pairs;
    for (const auto& wc_pair : wedge_count) {
        local_wedge_count_pairs.push_back(wc_pair);
    }
    
    // Serialize wedge count pairs
    vector<int> serialized_wedge_counts;
    for (const auto& pair : local_wedge_count_pairs) {
        serialized_wedge_counts.push_back(pair.first.first);
        serialized_wedge_counts.push_back(pair.first.second);
        serialized_wedge_counts.push_back(pair.second);
    }
    
    // Gather all serialized wedge counts
    int local_count_size = serialized_wedge_counts.size();
    vector<int> all_count_sizes(world_size);
    MPI_Allgather(&local_count_size, 1, MPI_INT, all_count_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    int total_count_size = 0;
    for (int size : all_count_sizes) {
        total_count_size += size;
    }
    
    vector<int> all_serialized_counts(total_count_size);
    vector<int> displacements(world_size, 0);
    for (int i = 1; i < world_size; i++) {
        displacements[i] = displacements[i-1] + all_count_sizes[i-1];
    }
    
    MPI_Allgatherv(serialized_wedge_counts.data(), local_count_size, MPI_INT,
                  all_serialized_counts.data(), all_count_sizes.data(), displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Merge wedge counts
    unordered_map<pair<Vertex, Vertex>, int, pair_hash> global_wedge_count;
    for (size_t i = 0; i < all_serialized_counts.size(); i += 3) {
        Vertex u1 = all_serialized_counts[i];
        Vertex u2 = all_serialized_counts[i+1];
        int count = all_serialized_counts[i+2];
        global_wedge_count[make_pair(u1, u2)] += count;
    }
    
    // Calculate butterfly counts
    for (const auto& wc_pair : global_wedge_count) {
        const pair<Vertex, Vertex>& endpoints = wc_pair.first;
        int d = wc_pair.second;
        Vertex u1 = endpoints.first;
        Vertex u2 = endpoints.second;
        butterfly_count[u1] += d * (d - 1) / 2;
        butterfly_count[u2] += d * (d - 1) / 2;
    }
    
    // Calculate center counts for vertices in my partition
    unordered_map<Vertex, int> center_count;
    for (const auto& wedge : my_wedges) {
        const pair<Vertex, Vertex>& endpoints = get<0>(wedge);
        Vertex center = get<2>(wedge);
        
        // Only count for centers in my partition
        if (vertex_partition.at(center) == rank) {
            center_count[center] += global_wedge_count[endpoints] - 1;
        }
    }
    
    // Merge center counts
    vector<pair<Vertex, int>> local_center_pairs;
    for (const auto& cc_pair : center_count) {
        local_center_pairs.push_back(cc_pair);
    }
    
    // Serialize center count pairs
    vector<int> serialized_center_counts;
    for (const auto& pair : local_center_pairs) {
        serialized_center_counts.push_back(pair.first);
        serialized_center_counts.push_back(pair.second);
    }
    
    // Gather all serialized center counts
    int local_center_size = serialized_center_counts.size();
    vector<int> all_center_sizes(world_size);
    MPI_Allgather(&local_center_size, 1, MPI_INT, all_center_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    int total_center_size = 0;
    for (int size : all_center_sizes) {
        total_center_size += size;
    }
    
    vector<int> all_serialized_centers(total_center_size);
    displacements.assign(world_size, 0);
    for (int i = 1; i < world_size; i++) {
        displacements[i] = displacements[i-1] + all_center_sizes[i-1];
    }
    
    MPI_Allgatherv(serialized_center_counts.data(), local_center_size, MPI_INT,
                  all_serialized_centers.data(), all_center_sizes.data(), displacements.data(),
                  MPI_INT, MPI_COMM_WORLD);
    
    // Add center counts to butterfly counts
    for (int i = 0; i < all_serialized_centers.size(); i += 2) {
        Vertex v = all_serialized_centers[i];
        int c = all_serialized_centers[i+1];
        butterfly_count[v] += c;
    }
    
    return butterfly_count;
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    string filename = "inputs/rMatGraph_J_5_500.txt";
    if (argc > 1) {
        filename = argv[1];
    }
    
    auto total_start = high_resolution_clock::now();
    
    // Only process 0 prints timing information
    auto load_start = high_resolution_clock::now();
    Graph G = loadLigraAdjacencyGraph(filename, world_rank, world_size);
    auto load_stop = high_resolution_clock::now();
    auto load_duration = duration_cast<milliseconds>(load_stop - load_start);
    
    auto preprocess_start = high_resolution_clock::now();
    auto G_prime = preprocess(G, [&](Vertex a, Vertex b) {
        return G.adj[a].size() > G.adj[b].size();
    });
    auto preprocess_stop = high_resolution_clock::now();
    auto preprocess_duration = duration_cast<milliseconds>(preprocess_stop - preprocess_start);

// METIS Partitioning
vector<int> part;
unordered_map<Vertex, int> vertex_partition;

if (world_rank == 0) {
    vector<Vertex> nodes = G_prime.vertices;
    int n = nodes.size();
    unordered_map<Vertex, idx_t> vertex_to_index;
    for (idx_t i = 0; i < n; ++i) {
        vertex_to_index[nodes[i]] = i;
    }

    // Build adjacency arrays for METIS
    vector<idx_t> xadj(n + 1, 0);
    vector<idx_t> adjncy;
    for (idx_t i = 0; i < n; ++i) {
        Vertex u = nodes[i];
        const auto& neighbors = G_prime.adj.at(u);
        for (Vertex v : neighbors) {
            adjncy.push_back(vertex_to_index[v]);
        }
        xadj[i + 1] = adjncy.size();
    }

    // METIS variables
    idx_t nvtxs = n;
    idx_t ncon = 1;
    idx_t nparts = world_size;
    idx_t objval;
    vector<idx_t> part_metis(n);
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

    // Partition the graph
    int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(),
                                  NULL, NULL, NULL, &nparts, NULL, NULL, options, &objval, part_metis.data());
    if (ret != METIS_OK) {
        cerr << "METIS partitioning failed." << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    part.assign(part_metis.begin(), part_metis.end());
}

// Broadcast partition vector
int n_part = (world_rank == 0) ? G_prime.vertices.size() : 0;
MPI_Bcast(&n_part, 1, MPI_INT, 0, MPI_COMM_WORLD);

if (world_rank != 0) {
    part.resize(n_part);
}
MPI_Bcast(part.data(), n_part, MPI_INT, 0, MPI_COMM_WORLD);

// Build vertex_partition map
vector<Vertex> nodes = G_prime.vertices;
for (int i = 0; i < nodes.size(); ++i) {
    vertex_partition[nodes[i]] = part[i];
}
    auto wedges_start = high_resolution_clock::now();
    auto wedges = getWedges(G_prime, world_rank, world_size, vertex_partition);
    auto wedges_stop = high_resolution_clock::now();
    auto wedges_duration = duration_cast<milliseconds>(wedges_stop - wedges_start);
    
    auto count_start = high_resolution_clock::now();
    auto butterfly_counts = countVertexButterflies(wedges, world_rank, world_size, vertex_partition);
    auto count_stop = high_resolution_clock::now();
    auto count_duration = duration_cast<milliseconds>(count_stop - count_start);
    
    auto total_stop = high_resolution_clock::now();
    auto total_duration = duration_cast<milliseconds>(total_stop - total_start);

    // Open output file for results
    ofstream ofs("mpi_results.txt");
    if (!ofs.is_open()) {
        cerr << "Error opening results file" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (world_rank == 0) {
        ofs << "\nVertex butterfly counts:\n";
        for (const auto& bc_pair : butterfly_counts) {
            ofs << "Vertex " << bc_pair.first << ": " << bc_pair.second << endl;
        }

        ofs << "\nTiming Information:\n";
        ofs << "--------------------------------\n";
        ofs << "Graph loading:      " << load_duration.count() << " ms\n";
        ofs << "Preprocessing:      " << preprocess_duration.count() << " ms\n";
        ofs << "Wedge computation:  " << wedges_duration.count() << " ms\n";
        ofs << "Butterfly counting: " << count_duration.count() << " ms\n";
        ofs << "--------------------------------\n";
        ofs << "Total time:         " << total_duration.count() << " ms\n";
        ofs << "Number of processes: " << world_size << "\n";
        ofs << "Total wedges: " << wedges.size() << endl;
        long long total_butterflies = 0;
        for (const auto& bc_pair : butterfly_counts) {
            total_butterflies += bc_pair.second;
        }
        ofs << "Total butterflies (sum of counts): " << total_butterflies << endl;
    }
    
    // Finalize MPI
    MPI_Finalize();
    return 0;
}
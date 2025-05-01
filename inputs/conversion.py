def convert_edge_list_to_adjacency_graph(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('%') and line.strip()]
    
    edges = []
    u_set, v_set = set(), set()

    for line in lines:
        u_str, v_str = line.strip().split()
        u = int(u_str)
        v = int(v_str)
        edges.append((u, v))
        u_set.add(u)
        v_set.add(v)

    # Normalize vertex IDs to start at 0
    u_list = sorted(list(u_set))
    v_list = sorted(list(v_set))
    u_map = {u: idx for idx, u in enumerate(u_list)}
    v_map = {v: idx for idx, v in enumerate(v_list)}
    num_u = len(u_list)
    num_v = len(v_list)

    # Build adjacency list for U side
    adj = [[] for _ in range(num_u)]
    for u, v in edges:
        adj[u_map[u]].append(v_map[v])

    # Build offset and edge list
    offsets = []
    edge_list = []
    count = 0
    for nbrs in adj:
        offsets.append(count)
        edge_list.extend(sorted(nbrs))
        count += len(nbrs)

    with open(output_file, 'w') as out:
        out.write("AdjacencyGraph\n")
        out.write(f"{num_u}\n")
        out.write(f"{len(edge_list)}\n")
        for offset in offsets:
            out.write(f"{offset}\n")
        for v in edge_list:
            out.write(f"{v}\n")

    print(f"✅ Converted to Ligra AdjacencyGraph format and saved to {output_file}")

# Example usage:
convert_edge_list_to_adjacency_graph("input", "my_input.txt")

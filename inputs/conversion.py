import random

def create_smaller_graph_file_by_edges(input_filename, output_filename, target_edge_count):
  """
  Reads a graph file, randomly selects a subset of edges, and writes it to a new file.

  Args:
      input_filename (str): The name of the input graph file.
      output_filename (str): The name of the output graph file.
      target_edge_count (int): The desired number of edges in the smaller graph.
  """

  edges = []
  with open(input_filename, 'r') as infile:
    infile.readline() # Skip the header
    for line in infile:
      line = line.strip()
      if not line:
        continue
      try:
        node1, node2 = map(int, line.split())
        edges.append((node1, node2))
      except ValueError:
        print(f"Skipping invalid line: {line}")

  # Randomly sample edges
  selected_edges = random.sample(edges, min(target_edge_count, len(edges))) # Ensure we don't sample more than we have

  with open(output_filename, 'w') as outfile:
    outfile.write("WeightedAdjacencyGraph\n")
    for node1, node2 in selected_edges:
      outfile.write(f"{node1} {node2}\n")

 #--- Main Execution ---
input_file = "rMatGraph_J_5_500.txt"
output_file = "smaller_graph_random.txt"

 # Calculate total edge count
with open(input_file, 'r') as infile:
    total_edge_count = len(infile.readlines()) - 1  # Subtract 1 for the header

target_edge_count = total_edge_count // 2  # Integer division to get approximately 1/4

create_smaller_graph_file_by_edges(input_file, output_file, target_edge_count)

print(f"Smaller graph file '{output_file}' created with approximately {target_edge_count} edges.")
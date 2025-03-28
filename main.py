import mysql.connector
import networkx
import community
import numpy
from sklearn.cluster import KMeans
import matplotlib.pyplot

conn = mysql.connector.connect(
    host="localhost",  # Usually "localhost" if running locally
    user="root",
    password="51050",
    database="gene_set_enrichment"
)
cursor = conn.cursor()

# Read from files to populate Genes, Pathways, and Gene_Pathway Tables
file_list = []
read_file = input("The file to read Pathway files from: ")
with open(read_file, "r") as f:
    for file_line in f:
          fparts = file_line.split()
          file_list.append(fparts[0])
          cursor.execute(
              """
              SELECT EXISTS(
                  SELECT 1 FROM Files 
                  WHERE fname = %s
              )
              """,
              (fparts[0],)
          )
          exists = cursor.fetchone()[0]
          # File name does not exist in Files table; insert into Files
          if not exists:
              if(len(fparts) > 1):
                cursor.execute(
                    "INSERT INTO Files (fname, scale) "
                    "VALUES (%s, %s)",
                    (fparts[0], fparts[1])
                )
              else:
                  cursor.execute(
                      "INSERT INTO Files (fname) "
                      "VALUES (%s)",
                      (fparts[0],)
                  )
          # File name exists in Files table
          else:
              if(len(fparts) > 1):
              # Update the new scale.
                  cursor.execute(
                    "UPDATE Files SET scale = %s WHERE fname = %s",
                    (fparts[1], fparts[0])
                  )
              # Else, nothing to update


# The 3-line code below analyzes ALL saved files
# We want to ONLY analyze listed files.
# ** ONE CONCERN : the MAX() normalizations will be incorrect if the MAX is in an older file not included in the current files to analyze **
# file_list = cursor.execute(
#   "SELECT fname FROM Files"
# ).fetchall()


for thisfile in file_list:
  with open(thisfile, "r") as f:
    for line in f:
        pathway, gene_list = line.strip().split("\t", 1)

        cursor.execute(
          '''
          SELECT EXISTS(
              SELECT 1
              FROM File_Pathway FP 
              JOIN Pathways P ON FP.pid = P.pid 
              WHERE FP.fname = %s AND P.pname = %s
          )
          ''',
          (thisfile, pathway)
        )
        exists = cursor.fetchone()[0]

        if not exists:
          # The pathway in this file does not exist
          # Populate Pathways Table (gid, pid, size)
          cursor.execute(
              "INSERT INTO Pathways (pname) "
              "VALUES (%s)",
              (pathway,)
          )
          pid = cursor.lastrowid
          # Populate File_Pathway Table (fname, pid)
          cursor.execute(
              "INSERT INTO File_Pathway (fname, pid) "
              "VALUES (%s, %s)",
              (thisfile, pid)
          )
        # Populate Genes (gid, gname)  *gname is unique.
        for gene in gene_list.split():
            # Populate Genes (gid, gname)  *gname is unique.
            cursor.execute(
                "INSERT IGNORE INTO Genes (gname) "
                "VALUES (%s)",
                (gene,)
            )
            # Populate Gene_Pathway (pid, gid)  *duplicate Pathway names are allowed.
            cursor.execute(
                "INSERT IGNORE INTO Gene_Pathway (pid, gid) "
                "SELECT P.pid, G.gid "
                "FROM Pathways P "
                "JOIN File_Pathway FP ON P.pid = FP.pid "
                "JOIN Genes G ON G.gname = %s "
                "WHERE P.pname = %s AND FP.fname = %s",
                (gene, pathway, thisfile)
            )
            conn.commit()

# It is possible to have genes that occur a very low number of times, which will skew the entire dataset. Filter out from Gene_Pathway, which is used to create the correlation table

# Create temporary table for genes to delete
cursor.execute(
    '''
    DROP TABLE IF EXISTS Gene_Filter
    '''
)
threshold = 20
cursor.execute(
    '''
    CREATE TABLE Genes_Filter AS
        SELECT gid
        FROM Gene_Pathway GP
        GROUP BY gid
        HAVING COUNT(*) <= %s   
    ''',
    (threshold,)
)
cursor.execute(
    '''
    DELETE FROM Gene_Pathway
    WHERE gid IN (
        SELECT gid
        FROM Genes_Filter
    )
    '''
)
cursor.execute(
    '''
    DROP TABLE Genes_Filter
    '''
)
# Assign pathway sizes based on the filtered Gene_Pathway Table
cursor.execute(
    '''
    UPDATE Pathways P 
    LEFT JOIN (
        SELECT pid, COUNT(*) AS psize
        FROM Gene_Pathway
        GROUP BY pid
    ) PS ON P.pid = PS.pid
    SET P.size = PS.psize
    '''
)
# Delete empty pathways that have been completely filtered out
cursor.execute(
    '''
    DELETE FROM Pathways
    WHERE size = 0 OR size IS NULL
    '''
)

# Get the max Files scale to be used in normalization
# If file_max_scale is 1.0, it means all scales are equal to 1.0
cursor.execute(
    '''
    SELECT MAX(scale) FROM Files
    '''
)
file_max_scale = cursor.fetchone()[0]
# Get the max pathway size to be used in normalization
cursor.execute(
    '''
    SELECT MAX(size) FROM Pathways
    '''
)
max_size = cursor.fetchone()[0] or 0

# Compute normalized Pathway weights
if max_size > 0:
    cursor.execute(
        '''
        UPDATE Pathways P
        JOIN File_Pathway FP ON P.pid = FP.pid
        JOIN Files F ON FP.fname = F.fname
        SET P.weight = (LOG10(%s + 1) - LOG10(P.size + 1)) * F.scale / %s
        ''',
        (max_size, file_max_scale)
    )
else:
    print("All pathways are empty.\n")
    exit(1)

# Create temporary table for calculating total occurrence weights for efficiency
cursor.execute(
    '''
    DROP TABLE IF EXISTS Gene_Total_Weight
    '''
)
cursor.execute(
    '''
    CREATE TABLE Gene_Total_Weight(
        SELECT GP.gid, SUM(P.weight) AS total_weight
        FROM Gene_Pathway GP
        JOIN Pathways P ON GP.pid = P.pid
        GROUP BY GP.gid
    )
    '''
)
# Create correlation matrix
# Filters correlation matrix: no entries with jaccard similarity below 0.3
cursor.execute(
    '''
    DROP VIEW IF EXISTS Correlation_Table
    '''
)
cursor.execute(
    '''
    CREATE VIEW Correlation_Table AS
    SELECT gp1.gid AS gene1, gp2.gid AS gene2,
    SUM(P.weight) / (GW1.total_weight + GW2.total_weight - SUM(P.weight)) AS jaccard_similarity
    FROM Gene_Pathway gp1
    JOIN Gene_Pathway gp2 ON gp1.pid = gp2.pid
    JOIN Pathways P ON P.pid = gp1.pid
    JOIN Gene_Total_Weight GW1 ON gp1.gid = GW1.gid
    JOIN Gene_Total_Weight GW2 ON gp2.gid = GW2.gid
    WHERE gp1.gid < gp2.gid
    GROUP BY gp1.gid, gp2.gid, GW1.total_weight, GW2.total_weight
    HAVING jaccard_similarity > 0.0
    '''
)

# Convert correlation matrix into network graph
cursor.execute(
    '''
    SELECT * FROM Correlation_Table
    '''
)
edges = cursor.fetchall()
G = networkx.Graph()
for gene1, gene2, similarity in edges:
    G.add_edge(gene1, gene2, weight = similarity)
partition = community.best_partition(G, weight = "weight")
# partition is a dict of genes in communities, ex) {Gene1: Community1}

with open("testing.txt", "w") as f:
    print("Graph nodes: ", G.nodes(), file=f)
    print("Graph edges: ", G.edges(data=True), file=f)

# Create a Louvain modularity community dict, ex) {Community1: Gene1, Gene2, ..}
louvain_communities = {}
for node, community_id in partition.items():
    if community_id not in louvain_communities:
        louvain_communities[community_id] = []
    louvain_communities[community_id].append(node)

final_clusters = {}

# Apply spectral clustering on these communities
for community_id, nodes in louvain_communities.items():
    subgraph = G.subgraph(nodes) # Subgraph from the community

    # Skip communities with less than 3 nodes (irrelevant clustering)
    if len(subgraph) < 3:
        for node in subgraph.nodes():
            final_clusters[node] = community_id
        continue

    # Compute the Laplacian Matrix (D - A) from the subgraph
    L = networkx.laplacian_matrix(subgraph, weight = "weight").toarray()
    # Compute eigenvectors and eigenvalues from the Laplacian Matrix
    eigenvalues, eigenvectors = numpy.linalg.eigh(L)

    # Initialize the best partition as the subgraph (the original graph in each Louvain community), the best number of clusters as 1, and the best modularity as the subgraph's initial modularity
    best_k = 1
    best_partition = {node: partition[node] for node in subgraph.nodes()}
    best_modularity = community.modularity(partition, subgraph, weight = "weight")

    # Greedy algorithm: test incremental numbers of clusters by KMeans, proceed if the new clusters have a higher modularity than the current
    for k in range(2, min(6, len(subgraph))):
        # 6 is empirically the ideal max number of spectral clusters
        X = eigenvectors[:, 1:k+1] # 1 ~ k eigenvectors (skipping 0)
        kmeans = KMeans(n_clusters = k, random_state = 42, n_init = 10).fit(X)
        # Create new partition dict based on new Kmean clusters
        # ex) Node: 3.1 (Community 3, Label 1)
        new_partition = {node: f"{community_id}.{label}" for node, label in zip(subgraph.nodes(), kmeans.labels_)}
        # New modularity of the new partition
        new_modularity = community.modularity(new_partition, subgraph, weight = "weight")

        if new_modularity >= best_modularity:
            best_k = k
            best_partition = new_partition
            best_modularity = new_modularity
    # Store final best partition
    final_clusters.update(best_partition)

print("Final clusters", final_clusters)

output_clusters = {}
for node, cluster_id in final_clusters.items():
    if cluster_id not in output_clusters:
        output_clusters[cluster_id] = []
    output_clusters[cluster_id].append(node)

for cluster, genes in output_clusters.items():
    print(f"Cluster {cluster}, Genes: {genes}")

# Make SQL commits
conn.commit()
conn.close()
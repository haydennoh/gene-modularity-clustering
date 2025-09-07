#!/usr/bin/env python3
import os
import getpass
from collections import defaultdict
import mysql.connector
import networkx as nx
import numpy as np
from sklearn.cluster import KMeans
import community


LOW_FREQ_THRESHOLD = int(os.getenv("LOW_FREQ_THRESHOLD", "20"))
JACCARD_THRESHOLD = float(os.getenv("JACCARD_THRESHOLD", "0.3"))
MAX_SPECTRAL_K = int(os.getenv("MAX_SPECTRAL_K", "6"))
SPECTRAL_NODE_CAP = int(os.getenv("SPECTRAL_NODE_CAP", "5000"))

MYSQL_HOST = os.getenv("MYSQL_HOST", "127.0.0.1")
MYSQL_USER = os.getenv("MYSQL_USER", "root")
MYSQL_PASSWORD = os.getenv("MYSQL_PASSWORD", "")
MYSQL_DB = os.getenv("MYSQL_DB", "gene_set_enrichment")

# UTILITY: Creates INDEX (B+ trees) to speed up queries
def ensure_index(conn, table_name, index_name, index_cols_sql):
    """
    Ensure an index_name exists on table_name.
    index_cols_sql: columns to create an index on
    """
    sql_check = """
        SELECT COUNT(*)
        FROM information_schema.statistics
        WHERE table_schema = DATABASE()
          AND table_name = %s
          AND index_name = %s
    """
    with conn.cursor() as cur:
        cur.execute(sql_check, (table_name, index_name))
        exists = cur.fetchone()[0]
        if not exists:
            cur.execute(f"CREATE INDEX {index_name} ON {table_name} {index_cols_sql}")

def create_helpful_indexes(conn):
    """
    Add PKs/unique keys/indexes that speed up frequent queries.
    """
    stmts = [
        # Files
        "ALTER TABLE Files ADD PRIMARY KEY (fname)",
        # Pathways
        "ALTER TABLE Pathways ADD UNIQUE KEY uq_pname (pname)",
        "ALTER TABLE Pathways ADD INDEX idx_size (size)",
        # Genes
        "ALTER TABLE Genes ADD UNIQUE KEY uq_gname (gname)",
        # File_Pathway
        "ALTER TABLE File_Pathway ADD PRIMARY KEY (fname, pid)",
        "ALTER TABLE File_Pathway ADD INDEX idx_fp_pid (pid)",
        "ALTER TABLE File_Pathway ADD INDEX idx_fp_fname (fname)",
        # Gene_Pathway
        "ALTER TABLE Gene_Pathway ADD PRIMARY KEY (pid, gid)",
        "ALTER TABLE Gene_Pathway ADD INDEX idx_gp_gid (gid)",
        "ALTER TABLE Gene_Pathway ADD INDEX idx_gp_pid (pid)",
    ]
    with conn.cursor() as cur:
        for sql in stmts:
            try:
                cur.execute(sql)
            except mysql.connector.Error:
                # already exists
                conn.rollback()
    conn.commit()

# INPUT PARSING HELPERS
def upsert_files(conn, file_spec_lines):
    """file_spec_lines: example -> pathway_file.txt [scale]"""
    files_to_upsert = []
    for line in file_spec_lines:
        parts = line.strip().split()
        if not parts:
            continue
        fname = parts[0]
        scale = parts[1] if len(parts) > 1 else None
        files_to_upsert.append((fname, scale))

    if files_to_upsert:
        with conn.cursor() as cur:
            cur.executemany(
                """
                INSERT INTO Files (fname, scale)
                VALUES (%s, COALESCE(%s, 1))
                ON DUPLICATE KEY UPDATE scale = VALUES(scale)
                """,
                files_to_upsert
            )
    return [fname for fname, _ in files_to_upsert]

def parse_pathway_files(file_list):
    """
    Returns:
      file_pathway_names: {fname: set(pname)}
      pathway_genes: {(fname, pname): [gene1, gene2, ...]}
    """
    file_pathway_names = defaultdict(set)
    pathway_genes = defaultdict(list)

    for fname in file_list:
        with open(fname, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                # Input format: "pathway<TAB>gene1 gene2 ...\n"
                parts = line.rstrip().split("\t", 1)
                if len(parts) == 1:
                    pname = parts[0]
                    genes_str = ""
                else:
                    pname, genes_str = parts
                file_pathway_names[fname].add(pname)
                if genes_str:
                    pathway_genes[(fname, pname)].extend(genes_str.split())
    return file_pathway_names, pathway_genes

def upsert_pathways_and_file_pathway(conn, file_pathway_names):
    """Upserts all unique pathway names."""
    all_pnames = {(pname,) for pset in file_pathway_names.values() for pname in pset}
    if all_pnames:
        with conn.cursor() as cur:
            cur.executemany(
                """
                INSERT INTO Pathways (pname)
                VALUES (%s)
                ON DUPLICATE KEY UPDATE pname = VALUES(pname)
                """,
                list(all_pnames)
            )
    # Build pid map
    with conn.cursor() as cur:
        cur.execute("SELECT pid, pname FROM Pathways")
        pid_by_pname = {pname: pid for pid, pname in cur.fetchall()}

    # Upsert File_Pathway pairs
    fp_pairs = []
    for fname, pset in file_pathway_names.items():
        for pname in pset:
            pid = pid_by_pname.get(pname)
            if pid is not None:
                fp_pairs.append((fname, pid))
    if fp_pairs:
        with conn.cursor() as cur:
            cur.executemany(
                """
                INSERT INTO File_Pathway (fname, pid)
                VALUES (%s, %s)
                ON DUPLICATE KEY UPDATE pid = VALUES(pid)
                """,
                fp_pairs
            )
    return pid_by_pname

def upsert_genes_and_gene_pathway(conn, pathway_genes, pid_by_pname):
    # Upsert unique genes
    all_genes = {g for glist in pathway_genes.values() for g in glist}
    if all_genes:
        with conn.cursor() as cur:
            cur.executemany(
                "INSERT IGNORE INTO Genes (gname) VALUES (%s)",
                [(g,) for g in all_genes]
            )

    # gid map
    with conn.cursor() as cur:
        cur.execute("SELECT gid, gname FROM Genes")
        gid_by_gname = {gname: gid for gid, gname in cur.fetchall()}

    # Bulk insert into Gene_Pathway
    gp_pairs = set()
    for (fname, pname), genes in pathway_genes.items():
        pid = pid_by_pname.get(pname)
        if pid is None:
            continue
        for g in genes:
            gid = gid_by_gname.get(g)
            if gid is not None:
                gp_pairs.add((pid, gid))
    if gp_pairs:
        with conn.cursor() as cur:
            cur.executemany(
                "INSERT IGNORE INTO Gene_Pathway (pid, gid) VALUES (%s, %s)",
                list(gp_pairs)
            )

# SQL BUILD
def build_current_files_temp(conn, file_list):
    with conn.cursor() as cur:
        cur.execute("DROP TEMPORARY TABLE IF EXISTS Current_Files")
        cur.execute("CREATE TEMPORARY TABLE Current_Files (fname VARCHAR(255) PRIMARY KEY)")
        cur.executemany(
            "INSERT INTO Current_Files (fname) VALUES (%s)",
            [(fn,) for fn in file_list]
        )

def filter_low_freq_genes(conn, threshold):
    with conn.cursor() as cur:
        cur.execute("DROP TEMPORARY TABLE IF EXISTS Genes_Filter")
        cur.execute(
            """
            CREATE TEMPORARY TABLE Genes_Filter AS
            SELECT gid
            FROM Gene_Pathway
            GROUP BY gid
            HAVING COUNT(*) <= %s
            """,
            (threshold,)
        )
        cur.execute(
            """
            DELETE GP FROM Gene_Pathway GP
            JOIN Genes_Filter GF USING (gid)
            """
        )
        cur.execute("DROP TEMPORARY TABLE IF EXISTS Genes_Filter")

def update_pathway_sizes_and_clean(conn):
    with conn.cursor() as cur:
        cur.execute(
            """
            UPDATE Pathways P
            LEFT JOIN (
                SELECT pid, COUNT(*) AS psize
                FROM Gene_Pathway
                GROUP BY pid
            ) PS ON P.pid = PS.pid
            SET P.size = PS.psize
            """
        )
        cur.execute("DELETE FROM Pathways WHERE size IS NULL OR size = 0")

def compute_norm_maxima(conn):
    """Compute maxima for files of current run."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT MAX(F.scale)
            FROM Files F
            JOIN Current_Files CF ON F.fname = CF.fname
            """
        )
        file_max_scale = cur.fetchone()[0] or 1
        cur.execute(
            """
            SELECT MAX(P.size)
            FROM Pathways P
            JOIN File_Pathway FP ON P.pid = FP.pid
            JOIN Current_Files CF ON FP.fname = CF.fname
            """
        )
        max_size = cur.fetchone()[0] or 0
    return file_max_scale, max_size

def update_weights_scoped(conn, file_max_scale, max_size):
    if max_size == 0:
        raise RuntimeError("All pathways are empty after filtering (max_size=0).")
    with conn.cursor() as cur:
        cur.execute(
            """
            UPDATE Pathways P
            JOIN File_Pathway FP ON P.pid = FP.pid
            JOIN Files F ON FP.fname = F.fname
            JOIN Current_Files CF ON F.fname = CF.fname
            SET P.weight = (LOG10(%s + 1) - LOG10(P.size + 1)) * F.scale / NULLIF(%s,0)
            """,
            (max_size, file_max_scale)
        )

def get_gene_total_weight(conn):
    with conn.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS Gene_Total_Weight")
        conn.commit()
        cur.execute(
            """
            CREATE TABLE Gene_Total_Weight AS
            SELECT GP.gid, SUM(P.weight) AS total_weight
            FROM Gene_Pathway GP
            JOIN Pathways P ON GP.pid = P.pid
            GROUP BY GP.gid
            """
        )
    ensure_index(conn, "Gene_Total_Weight", "idx_gtw_gid", "(gid)")

def get_correlation_table(conn, jaccard_threshold):
    with conn.cursor() as cur:
        cur.execute("DROP TABLE IF EXISTS Correlation_Table")
        conn.commit()
        cur.execute(
            """
            CREATE TABLE Correlation_Table AS
            SELECT gp1.gid AS gene1, gp2.gid AS gene2,
                   SUM(P.weight) / (GW1.total_weight + GW2.total_weight - SUM(P.weight)) AS jaccard_similarity
            FROM Gene_Pathway gp1
            JOIN Gene_Pathway gp2 ON gp1.pid = gp2.pid
            JOIN Pathways P ON P.pid = gp1.pid
            JOIN Gene_Total_Weight GW1 ON gp1.gid = GW1.gid
            JOIN Gene_Total_Weight GW2 ON gp2.gid = GW2.gid
            WHERE gp1.gid < gp2.gid
            GROUP BY gp1.gid, gp2.gid, GW1.total_weight, GW2.total_weight
            HAVING jaccard_similarity >= %s
            """,
            (jaccard_threshold,)
        )
    ensure_index(conn, "Correlation_Table", "idx_corr_g1", "(gene1)")
    ensure_index(conn, "Correlation_Table", "idx_corr_g2", "(gene2)")

def fetch_edges(conn):
    with conn.cursor() as cur:
        cur.execute("SELECT gene1, gene2, jaccard_similarity FROM Correlation_Table")
        return cur.fetchall()

# GRAPH CLUSTERING
def cluster_graph(edges):
    G = nx.Graph()
    for g1, g2, w in edges:
        G.add_edge(int(g1), int(g2), weight=float(w))

    if G.number_of_nodes() == 0:
        print("No edges above the Jaccard threshold; graph is empty.")
        return {}

    # Louvain macro-communities
    partition = community.best_partition(G, weight="weight")

    # Group by community id
    by_comm = defaultdict(list)
    for node, comm_id in partition.items():
        by_comm[comm_id].append(node)

    # Spectral clustering
    final_clusters = {}
    for comm_id, nodes in by_comm.items():
        subgraph = G.subgraph(nodes).copy()

        # Louvain label
        # Skip communities with less than 3 nodes (irrelevant clustering)
        if len(subgraph) < 3 or len(subgraph) > SPECTRAL_NODE_CAP:
            for n in subgraph.nodes():
                final_clusters[n] = str(comm_id)
            continue

        # Compute the Laplacian Matrix (D - A) from the subgraph
        L = nx.laplacian_matrix(subgraph, weight="weight").astype(float).toarray()
        eigenvalues, eigenvectors = np.linalg.eigh(L)

        # Start from Louvain labels
        best_partition = {n: str(comm_id) for n in subgraph.nodes()}
        best_modularity = community.modularity(best_partition, subgraph, weight="weight")

        # Try k=2..K
        for k in range(2, min(MAX_SPECTRAL_K, len(subgraph)) + 1):
            X = eigenvectors[:, 1:k]  # skip trivial eigenvector
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10).fit(X)
            new_partition = {
                n: f"{comm_id}.{label}"
                for n, label in zip(subgraph.nodes(), kmeans.labels_)
            }
            new_modularity = community.modularity(new_partition, subgraph, weight="weight")
            if new_modularity >= best_modularity:
                best_partition, best_modularity = new_partition, new_modularity

        final_clusters.update(best_partition)

    return final_clusters

def fetch_gid_to_name(conn):
  """Return {gid -> gname} map."""
  with conn.cursor() as cur:
      cur.execute("SELECT gid, gname FROM Genes")
      return {gid: gname for gid, gname in cur.fetchall()}

def persist_clusters_sql(conn, clusters, table_name="Gene_Communities"):
    """Create and fill a table storing (gid, label) for easy SQL querying."""
    with conn.cursor() as cur:
        cur.execute(f"DROP TABLE IF EXISTS {table_name}")
        cur.execute(
            f"""
            CREATE TABLE {table_name} (
                gid INT NOT NULL,
                label VARCHAR(64) NOT NULL,
                PRIMARY KEY (gid),
                FOREIGN KEY (gid) REFERENCES Genes(gid) ON DELETE CASCADE
            )
            """
        )
        rows = [(int(gid), str(label)) for gid, label in clusters.items()]
        if rows:
            cur.executemany(
                f"INSERT INTO {table_name} (gid, label) VALUES (%s, %s)",
                rows
            )

# MAIN
def main():
    # Build connection kwargs: prefer socket if provided
    conn_kwargs = {
        "host": MYSQL_HOST,
        "user": MYSQL_USER,
        "password": MYSQL_PASSWORD,
        "database": MYSQL_DB,
        "connection_timeout": 10,
    }

    # If no password and none in env, prompt interactively
    if not conn_kwargs.get("password") and os.getenv("MYSQL_PASSWORD") is None:
        conn_kwargs["password"] = getpass.getpass(f"MySQL password for {conn_kwargs.get('user','root')}: ")

    conn = mysql.connector.connect(**conn_kwargs)
    conn.autocommit = False

    try:
        create_helpful_indexes(conn)

        list_file = input("Path to file listing pathway files (and optional scales): ").strip()
        with open(list_file, "r") as f:
            file_spec_lines = [ln for ln in f if ln.strip()]

        # Files upsert
        file_list = upsert_files(conn, file_spec_lines)

        # Parse pathway files
        file_pathway_names, pathway_genes = parse_pathway_files(file_list)

        # Upsert Pathways & File_Pathway
        pid_by_pname = upsert_pathways_and_file_pathway(conn, file_pathway_names)

        # Upsert Genes & Gene_Pathway
        upsert_genes_and_gene_pathway(conn, pathway_genes, pid_by_pname)

        # Scope to current files
        build_current_files_temp(conn, file_list)

        # Filter low-frequency genes
        filter_low_freq_genes(conn, LOW_FREQ_THRESHOLD)

        # Update pathway sizes & clean empties
        update_pathway_sizes_and_clean(conn)

        # Normalize weights (scoped)
        file_max_scale, max_size = compute_norm_maxima(conn)
        if max_size == 0:
            raise RuntimeError("All pathways are empty after filtering.")
        update_weights_scoped(conn, file_max_scale, max_size)

        # get helper + correlation (with indexes)
        get_gene_total_weight(conn)
        get_correlation_table(conn, JACCARD_THRESHOLD)

        # Fetch edges & cluster
        edges = fetch_edges(conn)
        clusters = cluster_graph(edges)

        # Save clusters to output (manually written as output.txt)
        gid_to_name = fetch_gid_to_name(conn)

        # Group by cluster label
        grouped = defaultdict(list)
        for gid, label in clusters.items():
            grouped[label].append(gid)

        # Sort clusters in order of descending size
        items = sorted(grouped.items(), key=lambda kv: -len(kv[1]))

        with open("output.txt", "w") as out:
            total_assigned = sum(len(v) for v in grouped.values())
            out.write(f"Assigned genes: {total_assigned} (nodes in graph = {len(clusters)})\n\n")

            for label, gids in items:
                names = [gid_to_name.get(g, str(g)) for g in gids]
                out.write(f"Cluster {label} -> {len(names)} genes\n")
                out.write("  " + ", ".join(sorted(names)) + "\n\n")

        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()

if __name__ == "__main__":
    main()
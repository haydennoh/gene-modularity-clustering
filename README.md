# Gene Modularity Clustering

**gene-modularity-clustering** is a MySQL-integrated Python software library that produces statistically significant clusters of genes from existing pathways of diseases or cellular functions. The system constructs a correlation table from the database, generates **Louvain modularity communities** ("macro-clusters"), and applies **spectral clustering** for modularity-maximizing **"micro-clusters"**.

## Database Schema

The MySQL database schema consists of the following tables:

### `Genes`
| Column | Type | Constraints |
|--------|------|-------------|
| `gid`  | INT | PRIMARY KEY |
| `gname` | VARCHAR(255) | NOT NULL, UNIQUE |

### `Pathways`
| Column | Type | Constraints |
|---------|------|-------------|
| `pid`  | INT | PRIMARY KEY |
| `pname` | VARCHAR(255) | NOT NULL |
| `size` | INT | NOT NULL |
| `weight` | DOUBLE | |

### `Gene_Pathway`
| Column | Type | Constraints |
|---------|------|-------------|
| `pid`  | INT | PRIMARY KEY |
| `gid`  | INT | PRIMARY KEY |
| *PRIMARY KEY:* (`pid`, `gid`) |

### `Files`
| Column | Type | Constraints |
|---------|------|-------------|
| `fname` | VARCHAR(255) | PRIMARY KEY |
| `scale` | DOUBLE | Default: 1, must be `>= 1` |

### `File_Pathway`
| Column | Type | Constraints |
|---------|------|-------------|
| `fname` | VARCHAR(255) | PRIMARY KEY |
| `pid` | INT | PRIMARY KEY |
| *PRIMARY KEY:* (`fname`, `pid`) |

---

## Input Format

This library **requires** input files in a specific format. The user must provide a list of pathway files along with optional file scales.

### Example:
```
$ The input file to analyze: input.txt
```
#### `input.txt`
```
Datafile1.txt
Datafile2.txt 3
Datafile3.txt 4
Datafile4.txt
```
- `Datafile1.txt` (scale: default 1)
- `Datafile2.txt` (scale: 3)
- `Datafile3.txt` (scale: 4)
- `Datafile4.txt` (scale: default 1)

---

## Database Population & Normalization

### File Scale Normalization
The `Files` table is normalized based on the maximum scale in the dataset. Given the example where the **max scale** is **4**:
```
Datafile1.txt → scale = 1/4
Datafile2.txt → scale = 3/4
Datafile3.txt → scale = 1
Datafile4.txt → scale = 1/4
```

### Pathway Weight Normalization
The `Pathways` table weights are normalized based on the maximum pathway size. Given the maximum size is **9**:
```
Pathway1 → weight = 9/9 = 1
Pathway2 → weight = 8/9
Pathway3 → weight = 6/9
```
Each pathway’s weight is **multiplied by the file scale**:
```
Pathway1 (in Datafile1) → 1 * (1/4) = 1/4
Pathway2 (in Datafile1) → (8/9) * (1/4) = 2/9
Pathway3 (in Datafile1) → (6/9) * (1/4) = 3/18
```

---

## Correlation Matrix Computation

A **Jaccard similarity** correlation matrix is computed using weighted pathway occurrences:

\`
J(A, B) = |A ∩ B| / |A ∪ B|
\`

where:
- \( A \) is the sum of pathway weights that contain Gene A.
- \( B \) is the sum of pathway weights that contain Gene B.

This provides correlation values for unique **gene pairs**, normalized to prevent values exceeding **1**.

---

## Clustering Process

### **1. Louvain Modularity Community Formation**
- The **correlation matrix** is converted into a **weighted graph** where:
  - **Nodes** = Genes
  - **Edges** = Correlation values
- The **Louvain algorithm** partitions genes into **macro-communities**:
  ```
  Community1: Gene1, Gene2, Gene3
  Community2: Gene4, Gene5, ...
  ```
- These **macro-clusters** serve as the basis for further clustering.

### **2. Spectral Clustering for Micro-Clusters**
- Each **Louvain community** undergoes **spectral clustering**.
- The **Laplacian matrix (L)** is computed:
  \`
  L = D - A
  \`
  where:
  - \( D \) = Degree matrix (stores the number of edges per node)
  - \( A \) = Adjacency matrix (stores node connectivity)

Example:
```
D:         A:         L:
2  0  0    0  1  1    2 -1 -1
0  1  0    1  0  0   -1  1  0
0  0  1    1  0  0   -1  0  1
```

- **Eigenvector decomposition** identifies weakly connected subgroups.
- The **smallest nonzero eigenvalues** correspond to optimal splits.
- A **greedy algorithm** iteratively refines the clustering until modularity cannot be improved.

---

## Final Output

The algorithm produces **statistically significant gene clusters**, categorized into:
1. **Macro-clusters** (Louvain communities)
2. **Micro-clusters** (Spectral clustering refinements)

These clusters can be analyzed further using the **MySQL database**, retrieving additional details such as **common pathways, attributes, and relationships**.

---

## Installation

### **Dependencies**
- Python3
- MySQL
- Required Python libraries:
  ```bash
  pip install -r requirements.txt
  ```
  
### **Setup**
1. Clone the repository:
   ```bash
   git clone 
   cd gene-modularity-clustering
   ```
2. Configure MySQL database:
   ```sql
   CREATE DATABASE gene_cluster_db;
   ```
3. Run the setup script:
   ```bash
   python setup.py
   ```

---

## Usage

To analyze a dataset:
```bash
python main.py input.txt
```

Results are stored in the **MySQL database** and can be exported for further analysis.

---

## Author
- **Hayden Noh**

---


# Gene Modularity Clustering

**gene-modularity-clustering** is a Python + MySQL software library for constructing and clustering weighted gene networks from pathway data.  
It ingests pathway files, normalizes them by file scale and pathway size, computes a **weighted Jaccard similarity matrix** of genes, and applies a **two-stage clustering workflow**:  
1. **Louvain modularity** to detect macro-communities  
2. **Spectral K-means** refinement to extract micro-clusters with improved modularity  

Results are persisted in MySQL and exported as text summaries.

---

## Database Schema

The MySQL schema is optimized with **bulk upserts**, **B+-tree indexes**, and **scoped queries** for efficient scans.  

### `Genes`
| Column | Type | Constraints |
|--------|------|-------------|
| `gid`  | INT | PRIMARY KEY |
| `gname` | VARCHAR(255) | NOT NULL, UNIQUE |

### `Pathways`
| Column | Type | Constraints |
|---------|------|-------------|
| `pid`  | INT | PRIMARY KEY |
| `pname` | VARCHAR(255) | NOT NULL, UNIQUE |
| `size` | INT | |
| `weight` | DOUBLE | |

### `Gene_Pathway`
| Column | Type | Constraints |
|---------|------|-------------|
| `pid`  | INT | FOREIGN KEY → Pathways(pid) |
| `gid`  | INT | FOREIGN KEY → Genes(gid) |
| *PRIMARY KEY:* (`pid`, `gid`) |

### `Files`
| Column | Type | Constraints |
|---------|------|-------------|
| `fname` | VARCHAR(255) | PRIMARY KEY |
| `scale` | DOUBLE | Default: 1, must be `>= 1` |

### `File_Pathway`
| Column | Type | Constraints |
|---------|------|-------------|
| `fname` | VARCHAR(255) | FOREIGN KEY → Files(fname) |
| `pid` | INT | FOREIGN KEY → Pathways(pid) |
| *PRIMARY KEY:* (`fname`, `pid`) |

---

## Input Format

The library **requires** an input file listing pathway files and optional file-specific scales.

Example:

```bash
$ python main.py input.txt
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

#### Datafile1.txt
```
Pathway1    Gene1 Gene2 Gene3 Gene4 Gene5 Gene6
Pathway2    Gene4 Gene6 Gene9 Gene10 Gene11 Gene12 Gene13
Pathway3    Gene3 Gene7 Gene10 Gene13
```

---

## Database Population & Normalization

### File Scale Normalization
Each file’s scale is normalized by dividing by the **maximum scale among the files in the current run** (tracked in a temp table).  
Example with `max scale = 4`:

```
Datafile1.txt → scale = 1/4
Datafile2.txt → scale = 3/4
Datafile3.txt → scale = 1
Datafile4.txt → scale = 1/4
```

### Pathway Size Normalization
Pathway sizes are log-normalized relative to the largest pathway size in the current run:

$$
\
size(p) = \frac{\log_{10}(max size + 1)}{\log_{10}(|p| + 1)}
\
$$

This ensures smaller pathways carry proportionally higher weights.

### Pathway Weights
Final pathway weight:

$$
\
w(p) = size(p) \times file scale(f)
\
$$

### Gene Filtering
Genes that occur fewer than a threshold number of times (default: **20**, configurable via `LOW_FREQ_THRESHOLD`) are removed before building the correlation matrix.

---

## Correlation Matrix Computation

A **weighted Jaccard similarity** is computed between all gene pairs:

$$
\
J(A,B) = \frac{\sum_{p \in P(A) \cap P(B)} w(p)}{W(A) + W(B) - \sum_{p \in P(A) \cap P(B)} w(p)}
\
$$

where:
- \(P(A)\) = pathways containing gene A  
- \(W(A)\) = total pathway weight for gene A  

Only pairs with J(A,B) ≥ threshold (default **0.3**) are retained.

---

## Clustering Process

### 1. Louvain Modularity
- Graph built from correlation table:  
  - **Nodes** = Genes  
  - **Edges** = Jaccard similarity weights  
- Louvain algorithm partitions graph into **macro-communities**.  

### 2. Spectral Refinement
- For each Louvain community (with size ≥ 3 and ≤ cap, default **5000 nodes**):  
  - Build Laplacian matrix \( L = D - A \).  
  - Perform eigenvector decomposition.  
  - Apply K-means (k = 2…6 by default) on eigenvectors.  
  - Accept split only if modularity improves.  

---

## Final Output

The library outputs:  
1. **Macro-clusters** (Louvain)  
2. **Micro-clusters** (Spectral refinement)  

Results are:  
- Written to `output.txt` (gene lists grouped by cluster in order of largest size).  
- Stored in MySQL

---

## Installation

### Dependencies
- Python 3.9+  
- MySQL  
- Python libraries: `mysql-connector-python`, `networkx`, `numpy`, `scikit-learn`, `python-louvain`

### Setup
1. Clone repo:
   ```bash
   git clone https://github.com/your-repo/gene-modularity-clustering.git
   cd gene-modularity-clustering
   ```
2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Initialize MySQL schema:
   ```bash
   mysql -u root -p < setup.sql
   ```

---

## Usage

Run clustering:
```bash
python main.py input.txt
```

Configure parameters with environment variables (`LOW_FREQ_THRESHOLD`, `JACCARD_THRESHOLD`, `MAX_SPECTRAL_K`, etc.).

---

## Author
**Hayden Noh**

# Grid Prefix Filter
Grid Prefix Filter: Efficient Spatial Filtering for Approximate Nearest Neighbor Search

This repository contains the source code for the **Grid Prefix Filter** system, which performs approximate nearest neighbor (ANN) search for polygon datasets using a grid-based prefix filtering approach. The method applies geometric hashing and prefix sum techniques to efficiently prune irrelevant candidates, improving the scalability and accuracy of spatial similarity search based on Jaccard distance.

📌 **Overview**
Grid Prefix Filter uses a grid-based approach combined with a prefix filter to perform spatial filtering for polygon similarity search. The method first divides the space into a uniform grid, then filters out candidates based on their grid-based hash values, making the search more efficient.

- **Hashing Method**: Grid-based partitioning combined with prefix filtering.
- **Similarity Metric**: Jaccard distance based on geometric intersection area.
- **Pruning**: Efficient filtering through prefix sum techniques, reducing the number of candidates to process significantly.

📁 **Repository Structure**
```
GridPrefixFilter/
├── src/                  # C++ source code
│   ├── main.cpp          # Entry point
│   ├── mpi_filter.cpp    # MPI parallelism
│   ├── geoutil.cpp       # Geometry utilities
│   ├── query.cpp         # Prefix filter implementation
│   ├── brute_force.cpp   # Brute‐force implementation (exact method)
│   ├── util.h.cpp        # Basic utilities (other)
│   └── parse_geodata.cpp # Read WKT files
├── data/                 # Input polygon data (WKT format) and output CSVs (not included)
└── README.md             # This file
```

## ⚙️ **Dependencies**

* **GEOS (C API)** — tested with GEOS 3.12+
* **MPI (e.g., OpenMPI)**
* **C++20 compiler** (e.g., g++ or clang++)
* **Python** (optional for visualization)

  * `shapely`, `matplotlib` for Python-based analysis

---

## 🔧 **Building**

To build the Grid Prefix Filter system, run the following command:

```bash
make

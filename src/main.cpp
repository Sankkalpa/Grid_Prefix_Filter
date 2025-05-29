
#define GEOS_USE_ONLY_R_API
#include <unordered_set>  // implements std::unordered_set which stores elements in no particular order
#include <utility>
#include "parse_geodata.h"
#include "query.h"
#include "thread.h"
#include <cmath>
#include <list>
#include <mpi.h>
#include <iostream>
#include <sstream> // std::istringstream, std::ostringstream and std::stringstream for input and output operation on strings as if they were streams
#include <vector>  // implements std::vector to work on vectors
#include <map>     // implements std::map which stores key-value pairs in a sorted order based on keys 
#include <set>
#include <fstream> // std::ifstream for reading files, std::ofstream for writing of files and std::fstream for both
#include <string>  // implements std::string class for handling and manipulating 
#include <algorithm> // implements std::max, std::min, std::sort, std::find that can be applied to containers and sequences
#include <queue>   // implements std::queue, std::priority_queue
#include <numeric> // implements std::accumulate to sum up elements in a range 
#include <cstdlib> // std::exit and std::rand, for process control like EXIT_SUCCESS and EXIT_FAILURE 
#include <ctime>   // std::time and std::difftime to calculate differences in time

#include <geos_c.h>
#include "mpi_gis.h" 

//typedef std::unordered_map<std::string, QueryResult> AllResults;

using namespace std;

//Define some data types for readability
using ResultType = std::vector<std::pair<std::string, std::vector<std::pair<double, std::string>>>>;

using NeighborPair = std::pair<double, std::string>; // Represents distance and Neighbor ID

using QueryResult = std::vector<NeighborPair>; // All neighbors for a query


//prasing output from a file 
void ParseFromFile(GEOSContextHandle_t ctx, const std::string& filename, std::unordered_map<std::string, QueryResult>& allResults) {
    std::ifstream file(filename);
    std::string line;
    std::string currentQueryID;

    while (std::getline(file, line)) {
        // Check for Query ID line
        if (line.find("Query ID:") != std::string::npos) {
            // Extract the Query ID
            std::string idPrefix = "Query ID: ";
            size_t idStartPos = line.find(idPrefix) + idPrefix.length();
            currentQueryID = line.substr(idStartPos);
        }
        // Check for Neighbor line
        else if (!currentQueryID.empty() && line.find("Neighbor ID:") != std::string::npos) {
            std::string neighborPrefix = "Neighbor ID: ";
            size_t neighborStartPos = line.find(neighborPrefix) + neighborPrefix.length();
            size_t commaPos = line.find(", Distance:");
            std::string neighborID = line.substr(neighborStartPos, commaPos - neighborStartPos);

            std::string distancePrefix = "Distance: ";
            size_t distanceStartPos = line.find(distancePrefix) + distancePrefix.length();
            double distance = std::stod(line.substr(distanceStartPos));

            // Add the neighbor information to the current query's list
            allResults[currentQueryID].emplace_back(distance, neighborID);
        }
        // Reset currentQueryID if the line is empty, indicating a new query section
        else if (line.empty()) {
            currentQueryID.clear();
        }
    }

    file.close();
}

// -------------------------------------------------------------
// Helper structure for min-heap items (for top-k selection)
struct HeapItem {
    float sim;  // similarity value
    int index;  // training vector index
};

struct HeapComparator {
    bool operator()(const HeapItem& a, const HeapItem& b) {
        return a.sim > b.sim; // smallest similarity on top
    }
};

// -------------------------------------------------------------
// Compute Jaccard similarity for two vectors a and b.
template<typename T>
float jaccard_similarity(const vector<T>& a, const vector<T>& b) {
    float intersection = 0.0f, sum_a = 0.0f, sum_b = 0.0f;
    for (size_t i = 0; i < a.size(); i++) {
        float prod = a[i] * b[i];
        intersection += prod;
        sum_a += a[i];
        sum_b += b[i];
    }
    float union_val = sum_a + sum_b - intersection;
    return (union_val == 0.0f) ? 0.0f : (intersection / union_val);
}

// -------------------------------------------------------------
// Parse CSV files and process vectors.
// Each row: first value is an ID, remaining are vector components.
// Applies ceil-to-nearest-tenth and clips to [0,1].
void parse_vectors(int no_of_files,
                   const string& sketch_path,
                   vector<vector<float>>& all_vectors,
                   vector<int>& all_ids)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk_size = no_of_files / size;
    int remainder = no_of_files % size;
    int start = rank * chunk_size + min(rank, remainder);
    int end = (rank + 1) * chunk_size + (rank < remainder ? 1 : 0);

    if (end > no_of_files) end = no_of_files;
    // For each file (using a naming pattern)
    for (int i = start; i < end; i++) {
        string csv_file_path = sketch_path + to_string(i) + ".csv";
        ifstream csvfile(csv_file_path);
        if (!csvfile.is_open()) {
            cerr << "Could not open file: " << csv_file_path << endl;
            continue;
        }
        string line;
        while (getline(csvfile, line)) {
            stringstream ss(line);
            string cell;
            vector<string> cells;
            while (getline(ss, cell, ',')) {
                cells.push_back(cell);
            }
            if (cells.size() < 2) continue;
            int id = stoi(cells[0]);
            all_ids.push_back(id);
            vector<float> vec;
            for (size_t j = 1; j < cells.size(); j++) {
                float value = stof(cells[j]);
                // Ceil to nearest tenth
                value = ceil(value * 10.0f) / 10.0f;
                // Clip to [0,1]
                value = min(max(value, 0.0f), 1.0f);
                vec.push_back(value);
            }
            all_vectors.push_back(vec);
        }
        cout << csv_file_path << " has been read" << endl;
    }
    cout << "Finished processing vectors." << endl;
}

void parse_vectors_into_binary(int no_of_files,
                                 const string& sketch_path,
                                 vector<vector<int>>& all_vectors,
                                 vector<int>& all_ids)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunk_size = no_of_files / size;
    int remainder = no_of_files % size;
    int start = rank * chunk_size + min(rank, remainder);
    int end = (rank + 1) * chunk_size + (rank < remainder ? 1 : 0);
    // For each file (using a naming pattern)
    for (int i = start; i < end; i++) {
        // Construct the file name using the given naming pattern.
        string csv_file_path = sketch_path + to_string(i) + ".csv";
        ifstream csvfile(csv_file_path);
        if (!csvfile.is_open()) {
            cerr << "Could not open file: " << csv_file_path << endl;
            continue;  
        }
        
        string line;
        // int line_no = 1;
        while(getline(csvfile, line)) {
            stringstream ss(line);
            string cell;
            vector<string> cells;
            while(getline(ss, cell, ',')) {
                cells.push_back(cell);
            }
            if (cells.size() < 2)
                continue;
            
            int id = stoi(cells[0]);
            all_ids.push_back(id);
            
            vector<int> bin_vec;
            for (size_t j = 1; j < cells.size(); j++) {
                float value = stof(cells[j]);
                int bin_value = (value == 0.0f) ? 0 : 1;
                bin_vec.push_back(bin_value);
            }
            all_vectors.push_back(bin_vec);

            // cout << "Line number " << line_no <<" successfully read.\n";
            // line_no +=1;
        }
        cout << csv_file_path << " has been read." << endl;
    }
    cout << "Finished processing binary vectors." << endl;
}

//For query parsing

void parse_query(const string &query_file_path, vector<vector<float>> &vectors, vector<int> &ids) {
    ifstream csvfile(query_file_path);
    if (!csvfile.is_open()) {
        cerr << "Could not open query file: " << query_file_path << endl;
        return;
    }
    
    string line;
    while(getline(csvfile, line)) {
        stringstream ss(line);
        string cell;
        vector<string> cells;
        while(getline(ss, cell, ',')) {
            cells.push_back(cell);
        }
        if (cells.size() < 2) continue;
        
        // The first cell is the query id
        int id = stoi(cells[0]);
        ids.push_back(id);
        
        // Process the rest of the row into a vector
        vector<float> vec;
        for (size_t j = 1; j < cells.size(); j++) {
            float value = stof(cells[j]);
            // Apply ceiling to the nearest tenth
            value = ceil(value * 10.0f) / 10.0f;
            // Clip the value to [0, 1]
            value = min(max(value, 0.0f), 1.0f);
            vec.push_back(value);
        }
        vectors.push_back(vec);
    }
    
    cout << "Finished parsing query file: " << query_file_path 
         << ". Total queries parsed: " << vectors.size() << endl;
}

void parse_query_into_binary(const string &query_file_path, 
                             vector<vector<int>> &vectors, 
                             vector<int> &ids) {
    ifstream csvfile(query_file_path);
    if (!csvfile.is_open()) {
        cerr << "Could not open query file: " << query_file_path << endl;
        return;
    }
    
    string line;
    while(getline(csvfile, line)) {
        stringstream ss(line);
        string cell;
        vector<string> cells;
        
        // Split the line by commas.
        while(getline(ss, cell, ',')) {
            cells.push_back(cell);
        }
        
        // The line must have at least an ID and one value.
        if (cells.size() < 2)
            continue;
        
        // First cell is assumed to be the query ID.
        int id = stoi(cells[0]);
        ids.push_back(id);
        
        // Process each remaining cell and convert it into binary.
        vector<int> bin_vec;
        for (size_t j = 1; j < cells.size(); j++) {
            // Convert to float.
            float value = stof(cells[j]);
            // If the value is exactly 0, store 0; else, store 1.
            int bin_value = (value == 0.0f) ? 0 : 1;
            bin_vec.push_back(bin_value);
        }
        vectors.push_back(bin_vec);
    }
    
    cout << "Finished parsing query file into binary: " << query_file_path 
         << ". Total queries parsed: " << vectors.size() << endl;
}

// -------------------------------------------------------------
// Get the subset of query indices for the current MPI process.
// total: total number of queries, rank: current process rank, size: total processes.
vector<int> get_indices_for_process(int total, int rank, int size) {
    int base = total / size;
    int rem = total % size;
    int start = rank * base + min(rank, rem);
    int count = base + (rank < rem ? 1 : 0);
    vector<int> indices;
    for (int i = 0; i < count; i++) {
        indices.push_back(start + i);
    }
    return indices;
}

// -------------------------------------------------------------
// Brute-force kNN processing for a given set of query indices.
template<typename T>
map<int, vector<int>> process_chunk(const vector<int>& query_indices, 
                                      const vector<vector<T>>& query_vectors, 
                                      const vector<vector<T>>& training_vectors, 
                                      int k) {
    map<int, vector<int>> results;
    for (int idx : query_indices) {
        const vector<T>& q = query_vectors[idx];
        priority_queue<HeapItem, vector<HeapItem>, HeapComparator> heap;
        for (size_t j = 0; j < training_vectors.size(); j++) {
            float sim = jaccard_similarity(q, training_vectors[j]);
            if (heap.size() < static_cast<size_t>(k)) {
                heap.push({sim, static_cast<int>(j)});
            } else if (sim > heap.top().sim) {
                heap.pop();
                heap.push({sim, static_cast<int>(j)});
            }
        }
        // Extract sorted indices from the heap.
        vector<HeapItem> top_items;
        while (!heap.empty()) {
            top_items.push_back(heap.top());
            heap.pop();
        }
        sort(top_items.begin(), top_items.end(), [](const HeapItem& a, const HeapItem& b) {
            return a.sim > b.sim;
        });
        vector<int> top_indices;
        for (auto& item : top_items) {
            top_indices.push_back(item.index);
        }
        cout << "Brute-force: Query " << idx << " processed." << endl;
        results[idx] = top_indices;
    }
    return results;
}

// Process chunks with similarity values
// Instead of map<int, vector<int>>, return sim/id pairs:
template<typename T>
map<int, vector<pair<float,int>>> process_chunk_with_sims(
    const vector<int>&       query_indices,
    const vector<vector<T>>& query_vectors,
    const vector<vector<T>>& training_vectors,
    int                      k)
{
    map<int, vector<pair<float,int>>> results;
    for (int qid : query_indices) {
        const auto& q = query_vectors[qid];
        priority_queue<HeapItem, vector<HeapItem>, HeapComparator> heap;
        for (int j = 0; j < (int)training_vectors.size(); ++j) {
            float sim = jaccard_similarity(q, training_vectors[j]);
            if ((int)heap.size() < k) {
                heap.push({sim, j});
            } else if (sim > heap.top().sim) {
                heap.pop();
                heap.push({sim, j});
            }
        }
        // extract into a vector<pair<float,int>>
        vector<pair<float,int>> top;
        while (!heap.empty()) {
            auto it = heap.top(); heap.pop();
            top.emplace_back(it.sim, it.index);
        }
        // sort descending
        sort(top.begin(), top.end(),
             [](auto &a, auto &b){ return a.first > b.first; });
        
        cout << "Brute-force: Query " << qid << " processed." << endl;
        results[qid] = move(top);
    }
    return results;
}


// -------------------------------------------------------------
// Build an inverted index from the training vectors.
// For each dimension, store training indices where the value is nonzero.
template<typename T>
vector<vector<int>> create_inverted_list(const vector<vector<T>>& training_vectors) {
    int dim = training_vectors[0].size();
    vector<vector<int>> inv_list(dim);
    for (size_t idx = 0; idx < training_vectors.size(); idx++) {
        const vector<T>& vec = training_vectors[idx];
        for (int i = 0; i < dim; i++) {
            if (vec[i] != 0.0f) {
                inv_list[i].push_back(idx);
            }
        }
    }
    return inv_list;
}

// -------------------------------------------------------------
// Worker function for kNN search using the inverted list.
// For each query index in sub_query_indices, only training vectors that
// share a nonzero dimension (as indicated by the inverted index) are examined.
template<typename T>
pair<map<int, vector<int>>, vector<int>> worker_function(const vector<int>& sub_query_indices, 
                                                         const vector<vector<T>>& query_vectors, 
                                                         const vector<vector<T>>& training_vectors, 
                                                         const vector<vector<int>>& index, int k) {
    map<int, vector<int>> results;
    vector<int> seen_id_counts;
    for (int query_idx : sub_query_indices) {
        const vector<T>& q = query_vectors[query_idx];
        priority_queue<HeapItem, vector<HeapItem>, HeapComparator> heap;
        unordered_set<int> seen_ids;
        // Iterate over nonzero indices in the query vector.
        for (size_t i = 0; i < q.size(); i++) {
            if (q[i] != 0.0f) {
                for (int tid : index[i]) {
                    if (seen_ids.find(tid) == seen_ids.end()) {
                        seen_ids.insert(tid);
                        float sim = jaccard_similarity(q, training_vectors[tid]);
                        if (heap.size() < static_cast<size_t>(k)) {
                            heap.push({sim, tid});
                        } else if (sim > heap.top().sim) {
                            heap.pop();
                            heap.push({sim, tid});
                        }
                    }
                }
            }
        }
        vector<HeapItem> top_items;
        while (!heap.empty()) {
            top_items.push_back(heap.top());
            heap.pop();
        }
        sort(top_items.begin(), top_items.end(), [](const HeapItem& a, const HeapItem& b) {
            return a.sim > b.sim;
        });
        vector<int> top_indices;
        for (auto& item : top_items) {
            top_indices.push_back(item.index);
        }
        cout << "Inverted-index: Query " << query_idx << " processed." << endl;
        results[query_idx] = top_indices;
        seen_id_counts.push_back(seen_ids.size());
    }
    return make_pair(results, seen_id_counts);
}

// With sim versions of the code 

// You already have these:
// struct HeapItem { float sim; int index; };
// struct HeapComparator { bool operator()(const HeapItem &a, const HeapItem &b) { return a.sim > b.sim; } };
// template<typename T> float jaccard_similarity(const std::vector<T>&, const std::vector<T>&);

template<typename T>
std::pair<
  std::map<int,std::vector<std::pair<float,int>>>,
  std::vector<int>
>
worker_function_with_prefix_and_counts(
    const std::vector<int>&              query_indices,
    const std::vector<std::vector<T>>&   query_vectors,
    const std::vector<std::vector<T>>&   training_vectors,
    const std::vector<std::vector<int>>& inv_list,
    int                                   k,
    double                                t  // prefix threshold
) {
    std::map<int,std::vector<std::pair<float,int>>> results;
    std::vector<int> seen_counts;
    seen_counts.reserve(query_indices.size());

    for (int qidx : query_indices) {
        auto const &q = query_vectors[qidx];
        // 1) gather nonzero dims of q
        std::vector<int> qi;
        for (int d = 0; d < (int)q.size(); ++d)
            if (q[d] != T(0)) qi.push_back(d);
        int qlen = qi.size();

        // 2) prefix length p
        int p = qlen - (int)std::ceil(t * qlen) + 1;
        p = std::max(0, std::min(p, qlen));

        // 3) scan inverted-list on prefix dims
        std::priority_queue<HeapItem,std::vector<HeapItem>,HeapComparator> pq;
        std::unordered_set<int> seen_ids;
        for (int pi = 0; pi < p; ++pi) {
            int dim = qi[pi];
            for (int tid : inv_list[dim]) {
                if (!seen_ids.insert(tid).second) continue;
                float sim = jaccard_similarity(q, training_vectors[tid]);
                if ((int)pq.size() < k) {
                    pq.push({sim, tid});
                } else if (sim > pq.top().sim) {
                    pq.pop();
                    pq.push({sim, tid});
                }
            }
        }
        // record how many unique training vectors we examined
        seen_counts.push_back((int)seen_ids.size());

        // extract top-k
        std::vector<std::pair<float,int>> top;
        while (!pq.empty()) {
            top.emplace_back(pq.top().sim, pq.top().index);
            pq.pop();
        }
        std::reverse(top.begin(), top.end());
        results[qidx] = std::move(top);
    }

    return {std::move(results), std::move(seen_counts)};
}

// -------------------------------------------------------------
// Simple serialization: convert map<int, vector<int>> to a string.
string serialize_results(const map<int, vector<int>>& results) {
    ostringstream oss;
    for (auto &kv : results) {
        oss << kv.first;
        for (auto &v : kv.second)
            oss << " " << v;
        oss << "\n";
    }
    return oss.str();
}

// Deserialization: convert string back to map<int, vector<int>>
map<int, vector<int>> deserialize_results(const string& s) {
    map<int, vector<int>> res;
    istringstream iss(s);
    string line;
    while(getline(iss, line)) {
        istringstream ls(line);
        int key;
        ls >> key;
        vector<int> values;
        int val;
        while(ls >> val) {
            values.push_back(val);
        }
        res[key] = values;
    }
    return res;
}

//--------------------------------------------------------------
// Serialization and Deserialization of chunks for map data

// Serialize map<qID, vector<(sim,id)>> to a string
string serialize_sims(const map<int, vector<pair<float,int>>>& M) {
    ostringstream o;
    for (auto &kv : M) {
        int qid = kv.first;
        auto &vec = kv.second;
        o << qid << ' ' << vec.size();
        for (auto &pr : vec) {
            o << ' ' << pr.first << ' ' << pr.second;
        }
        o << '\n';
    }
    return o.str();
}

// Deserialize back
map<int, vector<pair<float,int>>> deserialize_sims(const string &s) {
    map<int, vector<pair<float,int>>> M;
    istringstream i(s);
    string line;
    while (getline(i, line)) {
        istringstream ls(line);
        int qid, cnt;
        ls >> qid >> cnt;
        vector<pair<float,int>> vec;
        for (int j = 0; j < cnt; ++j) {
            float sim; int id;
            ls >> sim >> id;
            vec.emplace_back(sim, id);
        }
        M[qid] = move(vec);
    }
    return M;
}



// // -------------------------------------------------------------
// // Main function using MPI for parallel processing.

// int main(int argc, char** argv) {
//     MPI_Init(&argc, &argv);
    
//     int rank, nprocs;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
//     if (argc != 4) {
//         if (rank == 0) {
//             std::cerr << "Incorrect number of arguments!\n";
//             std::cerr << "Usage: mpirun -n <num_processes> ./parallel_main "
//                          "<dataFile> <queryFile> <k>\n";
//         }
//         MPI_Finalize();
//         return 1;
//     }
    
//     std::string sketch_path  = argv[1];
//     std::string qsketch_path = argv[2];
//     int k = std::atoi(argv[3]);
    
//     // 1) Parse your training & query data (binary)
//     std::vector<std::vector<float>> training_vectors;
//     std::vector<int>               training_ids;


//     // parse_vectors_into_binary(100, sketch_path, training_vectors, training_ids);
//     parse_vectors(100, sketch_path, training_vectors, training_ids);
    
//     std::vector<std::vector<float>> query_vectors;
//     std::vector<int>               query_ids;
//     parse_query(qsketch_path, query_vectors, query_ids);
    
//     // 2) Build a full list of all query indices: 0,1,…,N-1
//     int Q = query_vectors.size();
//     std::vector<int> all_query_indices(Q);
//     std::iota(all_query_indices.begin(),
//               all_query_indices.end(),
//               0);
    
//     // -------------------------------------------------
//     // Brute‑Force with similarity scores
//     auto local_brute_sims = process_chunk_with_sims(
//         all_query_indices,
//         query_vectors,
//         training_vectors,
//         k
//     );
//     MPI_Barrier(MPI_COMM_WORLD);
//     if (rank == 0) std::cout << "All ranks finished brute‑force kNN.\n";

//     // Serialize & gather brute‑force
//     std::string local_brute_str = serialize_sims(local_brute_sims);
//     std::string all_brute_str;
//     if (rank == 0) {
//         all_brute_str = local_brute_str;
//         for (int p = 1; p < nprocs; ++p) {
//             int len;
//             MPI_Recv(&len, 1, MPI_INT, p, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             std::vector<char> buf(len+1);
//             MPI_Recv(buf.data(), len, MPI_CHAR, p, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             buf[len] = '\0';
//             all_brute_str += std::string(buf.data());
//         }
//     } else {
//         int len = local_brute_str.size();
//         MPI_Send(&len, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
//         MPI_Send(local_brute_str.data(), len, MPI_CHAR, 0, 100, MPI_COMM_WORLD);
//     }
    
//     // -------------------------------------------------
//     // Inverted‑Index with similarity scores
//     // build index
//     auto inv_list = create_inverted_list(training_vectors);
//     for (auto &v : inv_list)
//         if (v.empty()) v.push_back(0);
    
//     auto [local_index_sims, local_seen_counts] = worker_function_with_prefix_and_counts(
//         all_query_indices,
//         query_vectors,
//         training_vectors,
//         inv_list,
//         k,
//         0.8
//     );




//     MPI_Barrier(MPI_COMM_WORLD);
//     if (rank == 0) std::cout << "All ranks finished inverted‑index kNN.\n";

//     // Serialize & gather inverted‑index
//     std::string local_index_str = serialize_sims(local_index_sims);
//     std::string all_index_str;
//     if (rank == 0) {
//         all_index_str = local_index_str;
//         for (int p = 1; p < nprocs; ++p) {
//             int len;
//             MPI_Recv(&len, 1, MPI_INT, p, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             std::vector<char> buf(len+1);
//             MPI_Recv(buf.data(), len, MPI_CHAR, p, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             buf[len] = '\0';
//             all_index_str += std::string(buf.data());
//         }
//     } else {
//         int len = local_index_str.size();
//         MPI_Send(&len, 1, MPI_INT, 0, 200, MPI_COMM_WORLD);
//         MPI_Send(local_index_str.data(), len, MPI_CHAR, 0, 200, MPI_COMM_WORLD);
//     }
    
//     MPI_Barrier(MPI_COMM_WORLD);
//     if (rank == 0) std::cout << "Gathered all inverted‑index results.\n";

//     // -------------------------------------------------
//     // 4A) Deserialize & Merge Brute‑Force on rank 0
//     std::map<int, std::vector<std::pair<float,int>>> merged_brute, final_brute_topk;
//     if (rank == 0) {
//         std::istringstream brut_in(all_brute_str);
//         std::string line;
//         while (std::getline(brut_in, line)) {
//             auto part = deserialize_sims(line + "\n");
//             for (auto &kv : part) {
//                 auto &dst = merged_brute[kv.first];
//                 dst.insert(dst.end(), kv.second.begin(), kv.second.end());
//             }
//         }
//         // extract top‑k per query
//         for (auto &kv : merged_brute) {
//             int qid = kv.first;
//             auto &cands = kv.second;
//             std::priority_queue<HeapItem,
//                                 std::vector<HeapItem>,
//                                 HeapComparator> pq;
//             for (auto &pr : cands) {
//                 if ((int)pq.size() < k) {
//                     pq.push({pr.first, pr.second});
//                 } else if (pr.first > pq.top().sim) {
//                     pq.pop();
//                     pq.push({pr.first, pr.second});
//                 }
//             }
//             std::vector<std::pair<float,int>> top;
//             while (!pq.empty()) {
//                 top.emplace_back(pq.top().sim, pq.top().index);
//                 pq.pop();
//             }
//             std::reverse(top.begin(), top.end());
//             final_brute_topk[qid] = std::move(top);
//         }
//     }

//     // 4B) Deserialize & Merge Inverted‑Index on rank 0
//     std::map<int, std::vector<std::pair<float,int>>> merged_index, final_index_topk;
//     if (rank == 0) {
//         std::istringstream idx_in(all_index_str);
//         std::string line;
//         while (std::getline(idx_in, line)) {
//             auto part = deserialize_sims(line + "\n");
//             for (auto &kv : part) {
//                 auto &dst = merged_index[kv.first];
//                 dst.insert(dst.end(), kv.second.begin(), kv.second.end());
//             }
//         }
//         for (auto &kv : merged_index) {
//             int qid = kv.first;
//             auto &cands = kv.second;
//             std::priority_queue<HeapItem,
//                                 std::vector<HeapItem>,
//                                 HeapComparator> pq;
//             for (auto &pr : cands) {
//                 if ((int)pq.size() < k) {
//                     pq.push({pr.first, pr.second});
//                 } else if (pr.first > pq.top().sim) {
//                     pq.pop();
//                     pq.push({pr.first, pr.second});
//                 }
//             }
//             std::vector<std::pair<float,int>> top;
//             while (!pq.empty()) {
//                 top.emplace_back(pq.top().sim, pq.top().index);
//                 pq.pop();
//             }
//             std::reverse(top.begin(), top.end());
//             final_index_topk[qid] = std::move(top);
//         }
//     }

//     // -------------------------------------------------
//     // Recall on rank 0
//     if (rank == 0) {
//         int recall_count = 0, total_preds = 0;
//         for (auto &kv : final_index_topk) {
//             int qid = kv.first;
//             auto &preds = kv.second;
//             auto it = final_brute_topk.find(qid);
//             if (it == final_brute_topk.end()) continue;
//             std::unordered_set<int> truth;
//             for (auto &pr : it->second) truth.insert(pr.second);
//             for (auto &pr : preds) {
//                 ++total_preds;
//                 if (truth.count(pr.second)) ++recall_count;
//             }
//         }
//         std::cout << "Recall = "
//                   << recall_count << " / " << total_preds
//                   << " = " << (double)recall_count/total_preds
//                   << "\n";
//     }

//     // -------------------------------------------------------------
//     // --- 6) Average processed-vector count for prefix method ---
//     // local_seen_counts is length Q
//     std::vector<int> global_seen(Q,0);
//     MPI_Reduce(local_seen_counts.data(), global_seen.data(), Q,
//                MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

//     if (rank == 0) {
//         long long sum = 0;
//         for (int c: global_seen) sum += c;
//         double avg = double(sum) / Q;
//         std::cout << "[6] Avg vectors examined per query = "
//                   << avg << "\n";
//     }

//     MPI_Finalize();
//     return 0;
// }
int main(int argc, char **argv)
{
    int rank, size;
    int numSplit;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm = MPI_COMM_WORLD;
    char hostname[256];
    gethostname(hostname, 256);

    if (argc != 4)
    {
        if (rank == 0)
        {
            cerr << "Incorrect number of arguments!\n";
            cerr << "Usage: mpirun -n <num_processes> ./parallel_main <dataFile> <inputFile>\n";
        }
        MPI_Finalize();
        return 1;
    }
    string datafile = string(argv[1]);
    string inputfile = string(argv[2]);
    unsigned int numNeighbor = atoi(argv[3]);

    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    GEOSWKTWriter* wktWriter = GEOSWKTWriter_create_r(ctx);
    cout <<"error1" <<endl;
    vector<string> localFiles = splitFileMPI(datafile, comm);
    cout <<"error2" <<endl;
   // vector<GEOSGeometry *> *geos = read_wkt(localFiles[0], ctx); 
    unordered_map<string, GEOSGeometry*> *geos = read_wkt_with_id(localFiles[0], ctx); //need to store unique id as well
    cout <<"error3" <<endl;

    //print id + geometries
    if (geos != nullptr && !geos->empty())
    {
        for (const auto& pair : *geos) {
            const string& id = pair.first;
            GEOSGeometry* geometry = pair.second;

            cout << "ID: " << id << " - input Geometry: ";
            printCoordinates(ctx, geometry);
            cout << endl;
        }
    }
    else {
        cerr << "No geometries to print or failed to read geometries." << endl;
    }

    if (geos->empty())
    {
        cerr << "Error: No geometries to process in process " << rank << endl;
        MPI_Finalize();
        GEOS_finish_r(ctx);
        return 1;
    }

    //read query geometries and their unique id
    // vector<string> localFiles = splitFileMPI(datafile, comm);
    unordered_map<string, GEOSGeometry*> *querygeos = read_wkt_with_id(inputfile, ctx); 

    //printing querygeos
    // if (querygeos != nullptr && !querygeos->empty()){
    //     for (const auto& pair : *querygeos) {
    //         const string& id = pair.first;
    //         GEOSGeometry* geometry = pair.second;
        //
    //         cout << "ID: " << id << " - Query  Geometry: ";
    //         printCoordinates(ctx, geometry);
    //         cout << endl;
    //     }
    //     }
    //  else {
    //     cerr << "No geometries to print or failed to read geometries." << endl;
    // }

    // Vector to hold centered geometries
    unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID;

    // Iterate over the unordered_map, center each geometry, and store it with the identifier
    for (const auto& pair : *querygeos) {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            // Store the centered geometry with its identifier in the map
            centeredQueryGeosWithID[pair.first] = centeredGeo;
        } else {
            std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
        }
    }

    //reading querygeos
    //  for (const auto& item : centeredQueryGeosWithID) {
    //     std::cout << "ID: " << item.first << ", Centered Geometry: ";
    //     printCoordinates(ctx, item.second); // ctx is your GEOSContextHandle_t
    //     std::cout << std::endl;
    // }

    //perform lsh query
    int gridSize = 50;
    // int hashLength = 5;
    // ResultType LSHoutput;
    ResultType PrefixFilteroutput;

    //  LSHHashFunction(*geos, rank, centeredqueryGeos, numNeighbor, gridSize, hashLength, LSHoutput, ctx);

    //calling lsh function which performs lsh query process   
    //   LSHHashFunction(*geos, rank, centeredQueryGeosWithID, numNeighbor, gridSize, hashLength, LSHoutput, ctx);

    if (rank == 0)cout << "Executing PrefixFilterFunction..." << endl;
    
    // Perform the prefix filter function
    PrefixFilterFunction(*geos, rank, centeredQueryGeosWithID, numNeighbor, gridSize, PrefixFilteroutput, ctx);

    if (rank == 0)cout << "PrefixFilterFunction completed." << endl;

    MPI_Barrier(MPI_COMM_WORLD);

      //printing PrefixFilter output
    //    for (const auto& queryResult : PrefixFilteroutput) {
    //     const auto& queryInfo = queryResult.first; // Query ID and Geometry
    //     const auto& neighbors = queryResult.second; // Neighbors
    //     cout << "Query ID: " << queryInfo << endl;

    //     // Write each neighbor ID and distance
    //     for (const auto& neighborPair : neighbors) {
    //         const auto& distance = neighborPair.first;
    //         const auto& neighborInfo = neighborPair.second;
    //         cout << "Neighbor ID: " << neighborInfo << ", Distance: " << distance << endl;
    //     }
    // }

    //for writing output in a file
    std::string outputFilename = "Prefix_results_rank_" + std::to_string(rank) + ".txt";
    std::ofstream outFile("/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/results/Prefix_results_rank_" + std::to_string(rank) + ".txt", std::ios::out);

    if (!outFile.is_open()) {
        std::cerr << "Failed to open file for writing results." << std::endl;
        MPI_Finalize();
        return 1; // or handle the error appropriately
    }

    // Assuming 'outFile' is successfully opened as per your code

    for (const auto& queryResult : PrefixFilteroutput) {
        const std::string& queryID = queryResult.first;  // Query ID
        const auto& neighbors = queryResult.second;  // Vector of (distance, neighbor ID) pairs

        // Write the query ID
        outFile << "Query ID: " << queryID << std::endl;

        // Write each neighbor ID and distance
        for (const auto& neighbor : neighbors) {
            outFile << "Neighbor ID: " << neighbor.second << ", Distance: " << neighbor.first << std::endl;
        }

        outFile << std::endl; // Blank line for separating entries
    }

    outFile.close(); // Don't forget to close the file after writing

    //Barrier to synchronize all processes before reading files
    MPI_Barrier(MPI_COMM_WORLD);

    std::unordered_map<std::string, QueryResult> globalResults;
    if (rank == 0) {
        std::string baseDirectory = "/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/results/";
        std::string baseFilename = "Prefix_results_rank_";
        std::string fileExtension = ".txt";

        std::unordered_map<std::string, QueryResult> allResult;

        // Process each file individually
        for (int i = 0; i < size; ++i) {
            std::string filename = baseDirectory + baseFilename + std::to_string(i) + fileExtension;
            std::ifstream inFile(filename);

            if (inFile.is_open()) {
               // std::cout << "Reading contents from file: " << filename << std::endl;

                // Reset allResult for each file
                allResult.clear();

                // ParsefromFile function fills allResult with data from the current file
                ParseFromFile(ctx, filename, allResult);

                // Merge fileResults into globalResults
                // Merge fileResults into globalResults
                for (const auto& [queryID, neighbors] : allResult) {
                    // Insert new neighbors and existing ones into a temporary vector
                    std::vector<NeighborPair> combinedNeighbors = globalResults[queryID];
                    combinedNeighbors.insert(combinedNeighbors.end(), neighbors.begin(), neighbors.end());

                    // Sort combined neighbors by distance
                    std::sort(combinedNeighbors.begin(), combinedNeighbors.end(),
                              [](const NeighborPair& a, const NeighborPair& b) { return a.first < b.first; });

                    // Keep only the top 3 neighbors
                    if (combinedNeighbors.size() > numNeighbor) {
                        combinedNeighbors.resize(numNeighbor);
                    }

                    // Update globalResults with the top 3 neighbors
                    globalResults[queryID] = combinedNeighbors;
                }

                inFile.close();

            }  
            else {
                std::cerr << "Failed to open file: " << filename << std::endl;
            }
        }
    }
    // // After processing all files, print the combined global results
    // std::cout << "\nPrefix Global Results:\n";
    // for (const auto& [queryID, neighbors] : globalResults) {
    //     std::cout << "Prefix Query ID: " << queryID << std::endl;
    //     for (const auto& [distance, neighborID] : neighbors) {
    //         std::cout << "  Prefix Neighbor ID: " << neighborID << ", Distance: " << distance << std::endl;
    //     }
    //     std::cout << std::endl; // Separator for readability
    // }
    

    //uncomment from here

    //start brute fore search
    // vector<pair<GEOSGeometry*, vector<pair<double, GEOSGeometry*>>>> queryResults;
    vector<pair<string, vector<pair<double, string>>>> queryBFResults;
    BFquery(rank, centeredQueryGeosWithID, *geos, numNeighbor, queryBFResults);

    //  //printing BF query results
    //  for (const auto& queryResult : queryBFResults) {
    //     const auto& queryInfo = queryResult.first; // Query ID and Geometry
    //     const auto& neighbors = queryResult.second; // Neighbors
    //     cout << "BF Query ID: " << queryInfo << endl;
    //
    //     // Write each neighbor ID and distance
    //     for (const auto& neighbor : neighbors) {
    //         cout << " BF Neighbor ID: " << neighbor.second << ", Distance: " << neighbor.first << std::endl;
    //     }
    // }

    MPI_Barrier(MPI_COMM_WORLD);

    // Construct the output filename using the process rank
    std::string outputBFFilename = "BFquery_results_rank_" + std::to_string(rank) + ".txt";
    // Open an output file stream
    std::ofstream outBFFile("/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/results/BFquery_results_rank_" + std::to_string(rank) + ".txt", std::ios::out);

    if (!outBFFile.is_open())
    {
        std::cerr << "Failed to open file for writing results." << std::endl;
        MPI_Finalize();
        return 1; // or handle the error appropriately
    }

    for (const auto& resultPair : queryBFResults) {
        const string& queryID = resultPair.first;
        const auto& neighborPairs = resultPair.second;

        // Write query ID instead of geometry WKT
        outBFFile << "Query ID: " << queryID << std::endl;

        for (const auto& neighborPair : neighborPairs) {
            double distance = neighborPair.first;
            const string& neighborID = neighborPair.second;

            // Write neighbor ID instead of geometry WKT
            outBFFile << "Neighbor ID: " << neighborID << ", Distance: " << distance << std::endl;
        }
        outBFFile << std::endl;// Separator for readability; // Separator for readability
    }

    outBFFile.close(); 

    MPI_Barrier(MPI_COMM_WORLD);

    std::unordered_map<std::string, QueryResult> globalBFResults;
    if (rank == 0) {
        std::string baseDirectory = "/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/results/";
        std::string baseFilename = "BFquery_results_rank_";
        std::string fileExtension = ".txt";

        std::unordered_map<std::string, QueryResult> allBFResult;

        // Process each file individually
        for (int i = 0; i < size; ++i) {
            std::string filename = baseDirectory + baseFilename + std::to_string(i) + fileExtension;
            std::ifstream inFile(filename);

            if (inFile.is_open()) {
               // std::cout << "Reading contents from file: " << filename << std::endl;

                // Reset allResult for each file
                allBFResult.clear();

                // ParsefromFile function fills allResult with data from the current file
                ParseFromFile(ctx, filename, allBFResult);

                // Merge fileResults into globalResults
                // Merge fileResults into globalResults
                for (const auto& [queryID, neighbors] : allBFResult) {
                    // Insert new neighbors and existing ones into a temporary vector
                    std::vector<NeighborPair> combinedNeighbors = globalBFResults[queryID];
                    combinedNeighbors.insert(combinedNeighbors.end(), neighbors.begin(), neighbors.end());

                    // Sort combined neighbors by distance
                    std::sort(combinedNeighbors.begin(), combinedNeighbors.end(),
                              [](const NeighborPair& a, const NeighborPair& b) { return a.first < b.first; });

                    // Keep only the top 3 neighbors
                    if (combinedNeighbors.size() > numNeighbor) {
                        combinedNeighbors.resize(numNeighbor);
                    }

                    // Update globalResults with the top 3 neighbors
                    globalBFResults[queryID] = combinedNeighbors;
                }

                inFile.close();

            }  
            else {
                std::cerr << "Failed to open file: " << filename << std::endl;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // After processing all files, print the combined global results
    // std::cout << "\nCombined BF Global Results:\n";
    // for (const auto& queryResult : globalBFResults) {
    //     const auto& queryInfo = queryResult.first; // Query ID and Geometry
    //     const auto& neighbors = queryResult.second; // Neighbors
    //     std::cout << "BF Query ID: " << queryInfo << endl;
    
    //     // Write each neighbor ID and distance
    //     for (const auto& neighbor : neighbors) {
    //         std::cout << " BF Neighbor ID: " << neighbor.second << ", Distance: " << neighbor.first << std::endl;
    //     }
    // }

    // start comparing performance (currently used)
    if (rank == 0) {

        FILE* output = fopen("/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/comparisonResults.csv", "w");

        if (output == NULL) {
            cerr << "Failed to open file for writing.\n";
        } else {
            fprintf(output, "query_id,avg_recall,mse,false_positive_rate,total_lsh_neighbors,total_bf_neighbors,total_matched_neighbors\n");

            double totalOverallRecall = 0.0;
            double totalOverallBFRecall = 0.0;
            double totalMSE = 0.0;
            int totalCount = 0;
            int totalFalsePositives = 0;
            int totalBFNeighborCount = 0;  // Counter for all BF neighbors, including duplicates
            int totalLSHNeighborCount = 0; // Counter for all LSH neighbors, including duplicates
            int totalMatchedNeighborsCount = 0;

            for (const auto& [queryID, bfNeighbors] : globalBFResults) {

                // Check if the queryID exists in the LSH results
                // cout << "Hi" << endl;

                auto lshResultIt = globalResults.find(queryID);
                if (lshResultIt == globalResults.end()) continue;
                auto& lshNeighbors = lshResultIt->second;
                int correctNeighbors = 0;
                double mse = 0.0;
                int mseidx = 0;

                // Count all BF neighbors (including duplicates)
                totalBFNeighborCount += bfNeighbors.size(); // Add the number of neighbors for this query

                // Count all LSH neighbors (including duplicates)
                totalLSHNeighborCount += lshNeighbors.size(); // Add the number of neighbors for this query

                // cout << "Total BF neighbors: " << bfNeighbors.size() << endl;
                // cout << "Total Prefix neighbors: " << lshNeighbors.size() << endl;

                // Matching BF and LSH neighbors for comparison
                for (const auto& [distance, lshNeighborID] : lshNeighbors) {

                    // cout << "In the for loop" << endl;

                    bool isMatched = false;
                    for (const auto& [bfDistance, bfNeighborID] : bfNeighbors) {
                        if (lshNeighborID == bfNeighborID) {
                            // This neighbor is correct
                            correctNeighbors++;
                            isMatched = true;

                            // Calculate MSE for matching neighbors
                            mse += pow(bfDistance - distance, 2);
                            mseidx++;
                            break;
                        }
                    }

                    

                    // False positive if no match is found in BF results
                    if (!isMatched) {
                        totalFalsePositives++;
                    }
                }

                // cout << "Correct neighbors calculated: " << correctNeighbors << endl;
                

                if (!lshNeighbors.empty()) {
                    double queryRecall = static_cast<double>(correctNeighbors) / lshNeighbors.size();
                    totalOverallRecall += queryRecall;
                }

                if (!bfNeighbors.empty()) {
                    double queryBFRecall = static_cast<double>(correctNeighbors) / bfNeighbors.size();
                    totalOverallBFRecall += queryBFRecall;
                }

                if (mseidx > 0) {
                    mse /= mseidx;
                    totalMSE += mse;
                }

                totalCount++;
                totalMatchedNeighborsCount += correctNeighbors;
            }

            // Calculate averages
            double averageRecall = totalCount > 0 ? totalOverallRecall / totalCount : 0;
            double averageBFRecall = totalCount > 0 ? totalOverallBFRecall / totalCount : 0;
            double averageMSE = totalCount > 0 ? totalMSE / totalCount : 0;
            double falsePositiveRate = totalLSHNeighborCount > 0 ? static_cast<double>(totalFalsePositives) / totalLSHNeighborCount : 0;

            // Output the results
            cout << "Average Recall: " << averageRecall << endl;
            cout << "Average Recall by BF neighbors: " << averageBFRecall << endl;
            cout << "MSE: " << averageMSE << endl;
            cout << "Total Matched Neighbors: " << totalMatchedNeighborsCount << endl;
            cout << "Total Neighbors (including duplicates): " << totalLSHNeighborCount << endl;
            cout << "Total BF Neighbors (including duplicates): " << totalBFNeighborCount << endl;
            cout << "Total False Positives rate: " << falsePositiveRate << endl;

            // Additional console output...

            fclose(output);
        }
    }

    MPI_Finalize();
    GEOS_finish_r(ctx);
    return 0;
}



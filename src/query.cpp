// Implementation of the two query methods: LSH and brute force

#include "query.h"
#include <mpi.h>


double sketchSum(const std::vector<double>& sketch) {
    double sum = 0.0;
    for (double value : sketch) {
        sum += value;
    }
    return sum;
}

void queryThreadFunct(
                     unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
                     unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap,
                     vector<HashMap<std::string>*> hashMaps,
                     int numNeighbors, 
                     int hashLength, 
                     vector<unsigned int>& seeds, 
                     int rank,
                    // std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>>& LSHoutput,
                     ResultType& LSHoutput,
                     vector<vector<GEOSGeometry*>>& global_MBR_grid,
                      const vector<bool>& cellsToKeep) {

    // Create a new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    
   
    // Generate sketches for Query Geometries
   

unordered_map<string, pair<GEOSGeometry*, vector<double>>> QuerygeometrySketchMap;

// // Start Hashing time for polygon 
  double startTimeQuerysketch = MPI_Wtime();

// // // // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
// for (const auto& pair : centeredQueryGeosWithID) {
//     const string& id = pair.first; // Unique ID of the geometry
//     GEOSGeometry* geometry = pair.second; // Geometry pointer
// 
//     vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
//     if (sketchResult.empty()) {
//         cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
//     } else {
//         QuerygeometrySketchMap[id] = make_pair(geometry, sketchResult);
//     }
//   //  std::cout<< " - query Sketch Size: " << sketchResult.size() << std::endl;
//     
// }
// 
// // OVERALL time for querying 
//     double endTimeQuerySketch = MPI_Wtime();
//     double totalTimeQuerySketch = endTimeQuerySketch - startTimeQuerysketch;
//     
//  double maxTotalTimeQuerySketch;
// 
// // Use MPI_Reduce to find the maximum hashing time across all processes
// MPI_Reduce(&totalTimeQuerySketch, &maxTotalTimeQuerySketch, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
// 
// // The master process (rank 0) prints the maximum hashing time
// if (rank == 0) 
// {
//     std::cout << "Maximum time for query sketch creation across all processes: " << maxTotalTimeQuerySketch << " seconds" << std::endl;
// }
 


// Create and trim sketch for each geometry in the unordered_map centeredQueryGeosWithID
for (const auto& pair : centeredQueryGeosWithID) {
    const string& id = pair.first; // Unique ID of the geometry
    GEOSGeometry* geometry = pair.second; // Geometry pointer

    vector<double> fullSketch = sketch(ctx, &global_MBR_grid, geometry); // Generate the full sketch first
    vector<double> trimmedSketch;
    for (size_t i = 0; i < fullSketch.size(); ++i) {
        if (cellsToKeep[i]) {
            trimmedSketch.push_back(fullSketch[i]); 
            // Keep only the cells marked by cellsToKeep
        }
         
    }
   // std::cout << "ID: " << id << " - Trimmed Sketch Size: " << trimmedSketch.size() << std::endl;
    if (trimmedSketch.empty())
     {
        cerr << "Warning: Trimmed sketch generated for query geometry with ID " << id << " is empty." << endl;
    } else 
    {
        QuerygeometrySketchMap[id] = make_pair(geometry, trimmedSketch);
    }
}


// OVERALL time for querying 
    double endTimeQuerySketch = MPI_Wtime();
    double totalTimeQuerySketch = endTimeQuerySketch - startTimeQuerysketch;
    
 double maxTotalTimeQuerySketch;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeQuerySketch, &maxTotalTimeQuerySketch, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0) 
{
    std::cout << "Maximum time for query sketch creation across all processes: " << maxTotalTimeQuerySketch << " seconds" << std::endl;
}
 
// //printing sketch along with its id
//  for (const auto& item : QuerygeometrySketchMap) {
//      const string& id = item.first; // Unique ID
//      const auto& geometryAndSketch = item.second; // Pair of geometry pointer and sketch vector
//       const vector<double>& sketchResult = geometryAndSketch.second; // Extracting the sketch vector
//  // 
// //       cout << "ID: " << id << " - Sketch: ";
// //     for (double value : sketchResult) {
// //         cout << value << " ";
//     //}
//     cout << endl; // End of a single sketch
//  }
//  
//Similar process as before for generating sketches for query geometries

// Start Hashing time for polygon 
double startTimeQuery = MPI_Wtime();
double hashingStart, hashingEnd, hashingTotal = 0.0;
double queryingStart, queryingEnd, queryingTotal = 0.0;
double queryingSortEnd, queryingSortTotal =0.0;
double hashmap_key_search_start, hashmap_key_search_end, hashmap_key_search_total = 0.0;

size_t totalPolygons = 0;
//size_t totalHashes = 0;
for (const auto& item : QuerygeometrySketchMap)
 {
 
 const string& queryID = item.first;
  GEOSGeometry* queryGeo = item.second.first;
 const vector<double>& sketch = item.second.second;
// const vector<double>& querySketch = item.second.second; 

 //vector<pair<double, pair<string, GEOSGeometry*>>> nn; // Nearest neighbors
  vector<pair<double, string>> nn;  // Nearest neighbors storing only unique id's for preventing memory leak
 
 set<string> visited; // Set to keep track of visited IDs


for (unsigned int h = 0; h < hashMaps.size(); h++) {
    // size_t totalHashes = 0; 
      hashingStart = MPI_Wtime();
    //vector<int> hash = LSHHashs(const_cast<vector<double>*>(&sketch), hashLength, seeds[h]);
     std::vector<int> hash = LSHHash(const_cast<std::vector<double>*>(&sketch), hashLength, &seeds);
    

//Stop timing for hashing and accumulate the duration
        hashingEnd = MPI_Wtime();
        hashingTotal += hashingEnd - hashingStart;

        // Start timing for querying using MPI_Wtime
        queryingStart = MPI_Wtime();
 //  cout << "Hash vector for query hash map " << h << ": ";
//         for (int value : hash) {
//             cout << value << " ";
//         }
//         cout << endl;
//        
//         //counting the number of polygons associated with that lsh hash key
//         int count = hashMaps[h]->getValuesCount(hash);
//        // cout << "Number of polygons for hash in map " << ": " << count << endl; 
//        
//      //   Accumulate total polygons and hash count
//     totalPolygons += count;
//    // totalHashes++;
//      cout<<"total num of hashes" <<  totalPolygons << endl;
        
 hashmap_key_search_start = MPI_Wtime();
    const list<string>* bucket = hashMaps[h]->get(hash);
    hashmap_key_search_end = MPI_Wtime();
     hashmap_key_search_total += hashmap_key_search_end - hashmap_key_search_start;
    if (bucket != nullptr) {
   // cout << "bucket is not empty" <<endl;
        // Process elements in the bucket
        for (const string& neighborID : *bucket) {
            // Check if the geometry ID has been visited
            if (visited.find(neighborID) == visited.end()) {
                // Check if neighborID exists in geometrySketchMap
                auto neighborSketchIt = geometrySketchMap.find(neighborID);
                if (neighborSketchIt != geometrySketchMap.end()) {
                    visited.insert(neighborID); // Mark this ID as visited
                     GEOSGeometry* neighborGeo = neighborSketchIt->second.first;
                    double distance = jaccardDistance(ctx, queryGeo, neighborGeo);
                   // double distance = SketchJaccardDistance(ctx, queryGeo, neighborGeo);
                    
                    // cout<< "distance is "<<distance << endl;
                   // nn.push_back(make_pair(distance, make_pair(neighborID, neighborGeo))); 
                    nn.push_back(make_pair(distance, neighborID)); 
                  
                }
            }
        }
    }

// //for sketch similarity
// if (bucket != nullptr) {
//             for (const string& neighborID : *bucket) {
//                 if (visited.find(neighborID) == visited.end()) {
//                     auto neighborSketchIt = geometrySketchMap.find(neighborID);
//                     if (neighborSketchIt != geometrySketchMap.end()) {
//                         visited.insert(neighborID); // Mark this ID as visited
//                         const vector<double>& neighborSketch = neighborSketchIt->second.second;  // Using the sketch for the neighbor
//                         
//                         // Calculate sketch distance
//                         double sketchDistance = SketchJaccardDistance(sketch, neighborSketch);
// 
//                         // Store the sketch distance along with the neighbor ID
//                         nn.push_back(make_pair(sketchDistance, neighborID));
//                     }
//                 }
//             }
//         }
   // Stop timing for querying and accumulate the duration
        queryingEnd = MPI_Wtime();
        queryingTotal += queryingEnd - queryingStart;
}

        // Sort the nearest neighbors by distance
      //  sort(nn.begin(), nn.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
       // Sort the nearest neighbors by distance in descending order
        sort(nn.begin(), nn.end(), [](const pair<double, string>& a, const pair<double, string>& b) {
            return a.first > b.first; // Compare distances in descending order
        });


        // Keep only the required number of nearest neighbors
        if (nn.size() > numNeighbors) {
            nn.resize(numNeighbors);
        }
        

        // Store the nearest neighbors for the current query
        LSHoutput.push_back(make_pair(queryID, nn));
        
         // Stop timing for querying and accumulate the duration
        queryingSortEnd = MPI_Wtime();
        queryingSortTotal += queryingSortEnd - queryingStart;

}
//printing contents of lshoutput
for (const auto& result : LSHoutput) {
    const string& queryID = result.first;  // The query ID
    const auto& neighbors = result.second;  // The vector of nearest neighbors for this query

    std::cout << "Query ID: " << queryID << std::endl;  // Print the query ID

    // Check if there are any neighbors found
    if (neighbors.empty()) {
        std::cout << "  No neighbors found." << std::endl;
    } else {
        for (const auto& neighbor : neighbors) {
            double distance = neighbor.first;  // The distance to the neighbor
            const string& neighborID = neighbor.second;  // The neighbor's ID

            // Print the neighbor's ID and the distance
            std::cout << "  Neighbor ID: " << neighborID << ", Distance: " << distance << std::endl;
        }
    }

    std::cout << std::endl;  // Extra newline for better readability between query results
}

// // //for sketch similarity
// // // Initialize variables to store start and end times and total durations
// // double hashingStart, hashingEnd, hashingTotal = 0.0;
// // double queryingStart, queryingEnd, queryingTotal = 0.0;
// // 
// for (const auto& item : QuerygeometrySketchMap) {
//     const string& queryID = item.first;
//     const vector<double>& querySketch = item.second.second;  // Using the sketch for the query
//     
// 
//     vector<pair<double, string>> nn;  // Nearest neighbors storing only unique id's for preventing memory leak
//     set<string> visited; // Set to keep track of visited IDs
// 
//     for (unsigned int h = 0; h < hashMaps.size(); h++) {
//     
//     hashingStart = MPI_Wtime();
//     
//     vector<int> hash = LSHHashs(const_cast<vector<double>*>(&querySketch), hashLength, seeds[h]);
//         
//         // Print the hash values for debugging
// cout << "Hash values for query ID " << queryID << ": ";
// for(int val : hash) {
//     cout << val << " ";
// }
// cout << endl;
//         
//         // Stop timing for hashing and accumulate the duration
//         hashingEnd = MPI_Wtime();
//         hashingTotal += hashingEnd - hashingStart;
// 
//         // Start timing for querying using MPI_Wtime
//         queryingStart = MPI_Wtime();
//         const list<string>* bucket = hashMaps[h]->get(hash);
//         if (bucket != nullptr) {
//             for (const string& neighborID : *bucket) {
//                 if (visited.find(neighborID) == visited.end()) {
//                     auto neighborSketchIt = geometrySketchMap.find(neighborID);
//                     if (neighborSketchIt != geometrySketchMap.end()) {
//                         visited.insert(neighborID); // Mark this ID as visited
//                         const vector<double>& neighborSketch = neighborSketchIt->second.second;  // Using the sketch for the neighbor
//                         
//                         // Calculate sketch distance
//                         double sketchDistance = SketchJaccardDistance(querySketch, neighborSketch);
// 
//                         // Store the sketch distance along with the neighbor ID
//                         nn.push_back(make_pair(sketchDistance, neighborID));
//                     }
//                 }
//             }
//         }
//         // Stop timing for querying and accumulate the duration
//         queryingEnd = MPI_Wtime();
//         queryingTotal += queryingEnd - queryingStart;
//     }
// 
//     // Sort the nearest neighbors by distance
//     sort(nn.begin(), nn.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
// 
//     // Keep only the required number of nearest neighbors
//     if (nn.size() > numNeighbors) {
//         nn.resize(numNeighbors);
//     }
// 
//     // Store the nearest neighbors for the current query
//     LSHoutput.push_back(make_pair(queryID, nn));
// }

double maxhashingTotal;
double maxqueryingTotal;
double maxqueryingSortTotal;
double maxhashmap_key_search_total;
// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&hashingTotal, &maxhashingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
MPI_Reduce(&queryingTotal,&maxqueryingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
MPI_Reduce(&queryingSortTotal,&maxqueryingSortTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//MPI_Reduce(&hashmap_key_search_total,&maxhashmap_key_search_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);




if (rank == 0) 
{  // Only the process with rank 0 prints the timing
    std::cout << "Max LSH Hashing Time for query geometry: " << maxhashingTotal << " seconds" << std::endl;
   // std::cout << "Max time for particular lsh key search in hashmap for: " << hashmap_key_search_total << " seconds" << std::endl;
    std::cout << "Max Querying Time(only hash map search for neighbor) for query geometry: " << maxqueryingTotal << " seconds" << std::endl;
     std::cout << "Max Querying Time(only hash map search for neighbor) for query geometry after sorting: " << maxqueryingSortTotal << " seconds" << std::endl;
}


 // OVERALL time for querying 
    double endTimeQuery = MPI_Wtime();
    double totalTimeQuery = endTimeQuery - startTimeQuery;
    
 double maxTotalTimeQuery;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeQuery, &maxTotalTimeQuery, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0) 
{
    std::cout << "Maximum time for query(lsh hashing, querying hashmaps) across all processes: " << maxTotalTimeQuery << " seconds" << std::endl;
}


 
    // Clean up the GEOS context
    GEOS_finish_r(ctx);
    
 
}



//copy geometries + id's
vector<GEOSGeometry *> *copyGeometriesFromMap(const unordered_map<string, GEOSGeometry *> *original, GEOSContextHandle_t ctx) 
{
    vector<GEOSGeometry *> *copy = new vector<GEOSGeometry *>();
    for (const auto &pair : *original) {
        // Make a copy of the geometry
        GEOSGeometry *copyGeometry = GEOSGeom_clone_r(ctx, pair.second);
        copy->push_back(copyGeometry);
    }
    return copy;
}



// Function to print the contents of each HashMap in the vector
void printAllHashMaps(const std::vector<HashMap<std::string>*>& maps) {
    for (size_t i = 0; i < maps.size(); ++i) {
        std::cout << "Hash Map " << i << ":" << std::endl;
        // Call the print method on the current HashMap instance
        maps[i]->print(true); 
    }
}

// void printAllHashMaps(const std::vector<HashMap<std::string>*>& maps, std::ostream& os = std::cout) {
//     for (size_t i = 0; i < maps.size(); ++i) {
//         os << "Hash Map " << i << ":" << std::endl;
//         maps[i]->print(os, true);
//     }
// }

void PrefixFilterFunction(
    const unordered_map<string, GEOSGeometry*>& geos, // Changed to unordered_map
    int rank,
    const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
    unsigned int numNeighbor,
    int gridSize,
    ResultType& PrefixFilteroutput,
    GEOSContextHandle_t ctx
) {
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    
    //copy original geometries as when original geometries was used, it was resulting in segmentation fault
    unordered_map<string, GEOSGeometry*> geosCopy = copyGeometriesMap(geos, ctx);
  

   // Vector to hold centered geometries
    unordered_map<string, GEOSGeometry*> centeredGeosWithID;

  // Iterate over the unordered_map, center each geometry, and store it with its unique id
    for (const auto& pair : geos) {
    GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
    if (centeredGeo) {
        // Store the centered geometry with its identifier in the map
        centeredGeosWithID[pair.first] = centeredGeo;
    } else {
        std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
    }
   }

 vector<GEOSGeometry*> geometries;
for (const auto& pair : geosCopy) {
    GEOSGeometry* originalGeometry = pair.second;
    // Call the centerGeometry function on the original geometry
    GEOSGeometry* centeredGeometry = centerGeometry(ctx, originalGeometry);
    if (centeredGeometry != nullptr) {
        // Successfully centered the geometry, add it to the vector
        geometries.push_back(centeredGeometry);
         GEOSGeom_destroy_r(ctx, originalGeometry);
    } else {
        // Centering failed, handle the error accordingly
        cerr << "Warning: Failed to center geometry for ID " << pair.first << endl;
    }
}

 GEOSGeometry *global_MBR = GlobalMbr_parallel(&geometries, rank);
    
    
  //Creating grid over MBR in root only
  vector<vector<GEOSGeometry *>> global_MBR_grid = createGrid(ctx,
                                            global_MBR, gridSize); 



// Declare a map to store unique ID, geometry, and its sketch
unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap;
unordered_map<string, pair<GEOSGeometry*, vector<double>>> trimmedGeometrySketchMap;


//vector<vector<double>> localSketches;


//  //declare vector to store unique id, geometry and its sketch
 // unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap;
 // Start construction timing  for polygon MBR, grids, sketch 
  double startTimeSket = MPI_Wtime();
 // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
   for (const auto& pair : centeredGeosWithID) 
  {
    const string& id = pair.first; // Unique ID of the geometry
    GEOSGeometry* geometry = pair.second; // Geometry pointer

    vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
    if (sketchResult.empty()) 
    {
        cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
    } else 
    {
        geometrySketchMap[id] = make_pair(geometry, sketchResult);
    }
  }



// OVERALL time for construction 
    double endTimeSket = MPI_Wtime();
    double totalTimeSket = endTimeSket - startTimeSket;
    
double maxTotalTimeSket;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeSket, &maxTotalTimeSket, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0) {
    std::cout << "total time for sketch creation input after trimming: " << maxTotalTimeSket << " seconds" << std::endl;
}

// cout << "Error 1" << endl;
 
// for writing sketch vector output in a file
std::string outputSketchFilename = "sports_query1_sketch" + std::to_string(rank) + ".csv";  // Changed extension to .csv
std::ofstream outSKFile("/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/sketch/sports_query/" + outputSketchFilename, std::ios::out);

//to write sketch in a file in csv format
if (!outSKFile.is_open()) {
    std::cerr << "Failed to open file for writing results." << std::endl;
    MPI_Finalize();
    //return 1; // or handle the error appropriately
}

for (const auto& polygons : geometrySketchMap) {
    const std::string& polygonID = polygons.first;  // Query ID
    const auto& geometryAndketch = polygons.second;  // Pair of GEOSGeom_t* and vector of doubles

    // Write the query ID
    outSKFile << polygonID << ",";

    // Write the sketch values, comma-separated
   const std::vector<double>& sketch = geometryAndketch.second;
    for (size_t i = 0; i < sketch.size(); ++i) {
        outSKFile << sketch[i];
        if (i != sketch.size() - 1) {
            outSKFile << ",";  // Add a comma only if it's not the last element
        }
    }

    outSKFile << "\n"; // New line for separating entries (removed extra std::endl for blank line)
}

outSKFile.close(); // Don't forget to close the file after writing

// cout << "Error 2" << endl;

// Creating inverted index
unordered_map<int, std::vector<std::string>> invertedIndex; // Inverted index map: key -> grid cell index, value -> list of geometry IDs

// Reserve space for the inverted index;i++ )

ulong totalVectors = 0;

for (const auto& pair : geometrySketchMap) {
    const string& id = pair.first;
    const vector<double>& sketch = pair.second.second;
    for (size_t i = 0; i < sketch.size(); ++i) {
        if (sketch[i] > 0.0) {
        //     int gridX = i % gridSize;
        //     int gridY = i / gridSize;
        //     int flattenedIndex = gridY * gridSize + gridX;

        //     cout << "Flattened Index: " << flattenedIndex << endl;

            // if (flattenedIndex >= gridSize * gridSize) {
            //     std::cerr << "Error: Flattened index out of bounds for geometry ID " << id << std::endl;
            // }
            // else {
            invertedIndex[i].push_back(id);
            }
            
        // }
    }

    totalVectors++;
}

// std::cout << "Size of inverted index" << invertedIndex.size() << endl;

// cout << "Error 3" << endl;

// Store the length of invertedIndex for each grid cell
unordered_map<int, int> invertedIndexLength;
for (int i = 0; i < gridSize * gridSize; ++i) {
    // If the grid cell is not populated in invertedIndex, initialize the length to 0
    if (invertedIndex.find(i) == invertedIndex.end()) {
        invertedIndexLength[i] = 0;
    } else {
        invertedIndexLength[i] = invertedIndex[i].size();  // Store the number of geometries in this cell
    }
}

// cout << "Error 4" << endl;


MPI_Barrier(MPI_COMM_WORLD); // Synchronize all processes before printing

// cout << "Error 5" << endl;

// Save the inverted index lengths to a file
std::string outputFilename = "inverted_index_lengths" + std::to_string(rank) + ".csv";  // Changed extension to .csv
std::ofstream outFile("/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/inverted_index_lengths/" + outputFilename, std::ios::out);
if (!outFile.is_open()) {
    std::cerr << "Failed to open file for writing inverted index lengths." << std::endl;
    MPI_Finalize();
    //return 1; // or handle the error appropriately
}
outFile << "index, length\n"; // Header for the CSV file
// Write the inverted index lengths to the file
for (size_t i = 0; i < invertedIndexLength.size(); ++i) {
    outFile << i << "," << invertedIndexLength[i] << endl;
    // if (i != invertedIndexLength.size() - 1) {
    //     outFile << ",";  // Add a comma only if it's not the last element
    // }
}
outFile.close(); // Don't forget to close the file after writing

// cout << "Error 6" << endl;



// To store the length of the seen_ids set
std::vector<unsigned int> seen_id_Lengths;  // Set to store the IDs of already processed geometries

// Start the timing for the query processing
double startTimeQuery = MPI_Wtime();


// Now, let's process each query:

// cout << "Starting to process each query..." << endl;
    for (const auto& queryPair : centeredQueryGeosWithID) {
        const string& queryID = queryPair.first;  // Query geometry ID
        GEOSGeometry* queryGeometry = queryPair.second;  // Query geometry pointer

        // cout << "Error 20" << endl;

        // Generate sketch for the query geometry (trimmed if needed)
        vector<double> querySketch = sketch(ctx, &global_MBR_grid, queryGeometry);

        // cout << "Error 21" << endl;

        double querySketchSum = sketchSum(querySketch);
        double traversalSum = querySketchSum;

        // cout << "Error 22" << endl;

        // Result vector for storing the top-k neighbors (Jaccard distances and geometry IDs)
        std::vector<std::pair<double, std::string>> queryResults;

        // cout << "Error 23" << endl;

        // Map to store non-zero grid cells along with their corresponding flattened index
        std::vector<int> nonZeroCellIndexes;

        // cout << "Processing query geometry ID: " << queryID << endl;

        // Iterate over each non-zero grid cell in the query sketch
        for (size_t i = 0; i < querySketch.size(); ++i) {
            if (querySketch[i] > 0.0) {  // If the value is non-zero
                // Calculate the flattened index for this grid cell
                // int gridX = i % gridSize;
                // int gridY = i / gridSize;
                // int flattenedIndex = gridY * gridSize + gridX;
            
                // Add the flattened index to the list of non-zero cell indexes
                nonZeroCellIndexes.push_back(i);
            }
        }

        // cout << "Error 24" << endl;

        // Sort the non-zero cell indexes in ascending order of the inverted list size
        std::sort(nonZeroCellIndexes.begin(), nonZeroCellIndexes.end(),
            [&invertedIndex](int a, int b) {
                return invertedIndex[a].size() < invertedIndex[b].size(); // Compare based on list size
            });

        // cout << "Error 25" << endl;

        // Iterate over the sorted non-zero grid cell indexes
        // Inside your function, define the seen_ids set
        std::unordered_set<std::string> seen_ids; // Set to store the IDs of already processed geometries

        for (int flattenedIndex : nonZeroCellIndexes) {

            // if (seen_ids.size() > 0.6 * totalVectors) {
            //     break;  // Stop if we have seen enough IDs
            // }

            
            // Look up the inverted index to find geometries that have non-zero values in this cell
            if (invertedIndex.find(flattenedIndex) != invertedIndex.end()) {
                const std::vector<std::string>& relevantGeometries = invertedIndex[flattenedIndex];

                // cout << "Error 25.1" << endl;

                // For each geometry in this grid cell, calculate the Jaccard distance
                for (const auto& geoID : relevantGeometries) {
                    // Check if this geometry has already been processed
                    if (seen_ids.find(geoID) != seen_ids.end()) {
                        continue;  // Skip if this geometry has already been processed
                    }

                    // cout << "Error 25.2" << endl;

                    // Get the corresponding geometry's sketch from geometrySketchMap
                    const auto& geoPair = geometrySketchMap[geoID];

                    // cout << "Error 25.3" << endl;

                    GEOSGeometry* geometry = geoPair.first;

                    // cout << "Error 25.4" << endl;

                    // Compute Jaccard distance between the query and the geometry
                    double jaccardDist = jaccardDistance(ctx, queryGeometry, geometry);

                    // cout << "Error 25.5" << endl;

                    // Add the result (distance, geometryID) to the queryResults vector
                    queryResults.push_back(std::make_pair(jaccardDist, geoID));

                    // cout << "Error 25.6" << endl;

                    // Mark this geometry as processed by adding its ID to seen_ids
                    seen_ids.insert(geoID);

                    // if (seen_ids.size() > 0.6 * totalVectors) {
                    //     break;  // Stop if we have seen enough IDs
                    // }

                    // cout << "Error 25.7" << endl;
                }

        
                // cout << "Error 25.8" << endl;
        
            }

            // cout << "Error 26" << endl;

            traversalSum -= querySketch[flattenedIndex];  // Decrease the sum of the query sketch

            if (traversalSum <= 0.7 * querySketchSum) {
                break; 
            }


        }

        // Sort the results by Jaccard distance (ascending order)
        std::sort(queryResults.begin(), queryResults.end(),
                  [](const std::pair<double, std::string>& a, const std::pair<double, std::string>& b) {
                      return a.first < b.first;  // Sort by the distance (smaller distances come first)
                });

    

        // Store the size of the seen_ids set in the seen_id_Lengths vector
        seen_id_Lengths.push_back(seen_ids.size());

    
        // Keep only the top-k nearest neighbors
        if (queryResults.size() > numNeighbor) {
            queryResults.resize(numNeighbor);
        }

        // Store the query ID and the list of top-k results in PrefixFilteroutput
        PrefixFilteroutput.push_back(std::make_pair(queryID, queryResults));


    }


// ulong querySize = PrefixFilteroutput.size();
// ulong globalSeenIDLengths = 0;
// MPI_Reduce(&seen_id_Lengths, &globalSeenIDLengths, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

// if (rank == 0) {
//     double avgSeenIDLengths = static_cast<double>(globalSeenIDLengths) / querySize;
//     std::cout << "Average number of seen IDs per query: " << avgSeenIDLengths << std::endl;
// }

// cout << "Error 7" << endl;

// End the timing for the query processing
double endTimeQuery = MPI_Wtime();
double totalTimeQuery = endTimeQuery - startTimeQuery;
double maxTotalTimeQuery;

// // Use MPI_Reduce to find the maximum query processing time across all processes
MPI_Reduce(&totalTimeQuery, &maxTotalTimeQuery, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// // The master process (rank 0) prints the maximum query processing time
if (rank == 0) {
    std::cout << "Maximum Query Processing in Spatial Prefix Filter: " << maxTotalTimeQuery << " seconds" << std::endl;
}

// cout << "Error 8" << endl;

// Save the seen_id_Lengths to a file
std::string outputSeenIDLengthsFilename = "seen_id_lengths" + std::to_string(rank) + ".csv";  // Changed extension to .csv
std::ofstream outSeenIDLengthsFile("/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/seen_id_lengths/" + outputSeenIDLengthsFilename, std::ios::out);
if (!outSeenIDLengthsFile.is_open()) {
    std::cerr << "Failed to open file for writing seen ID lengths." << std::endl;
    MPI_Finalize();
    //return 1; // or handle the error appropriately
}
// Write the seen_id_Lengths to the file
for (size_t i = 0; i < seen_id_Lengths.size(); ++i) {
    outSeenIDLengthsFile << seen_id_Lengths[i];
    if (i != seen_id_Lengths.size() - 1) {
        outSeenIDLengthsFile << ",";  // Add a comma only if it's not the last element
    }
}
outSeenIDLengthsFile << "\n"; // New line for separating entries
outSeenIDLengthsFile.close(); // Don't forget to close the file after writing

// cout << "Error 9" << endl;

}




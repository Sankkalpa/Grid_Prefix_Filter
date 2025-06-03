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

void PrefixFilterFunction(
    const unordered_map<string, GEOSGeometry*>& geos,
    int rank,
    const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
    unsigned int numNeighbor,
    int gridSize,
    ResultType& PrefixFilteroutput,
    GEOSContextHandle_t ctx
) {
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    // Copy original geometries to avoid segmentation fault when using originals
    unordered_map<string, GEOSGeometry*> geosCopy = copyGeometriesMap(geos, ctx);

    // Vector to hold centered geometries
    unordered_map<string, GEOSGeometry*> centeredGeosWithID;
    // Iterate over the unordered_map, center each geometry, and store it with its unique ID
    for (const auto& pair : geos) {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            centeredGeosWithID[pair.first] = centeredGeo;
        } else {
            std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
        }
    }

    vector<GEOSGeometry*> geometries;
    for (const auto& pair : geosCopy) {
        GEOSGeometry* originalGeometry = pair.second;
        GEOSGeometry* centeredGeometry = centerGeometry(ctx, originalGeometry);
        if (centeredGeometry != nullptr) {
            geometries.push_back(centeredGeometry);
            GEOSGeom_destroy_r(ctx, originalGeometry);
        } else {
            cerr << "Warning: Failed to center geometry for ID " << pair.first << endl;
        }
    }

    GEOSGeometry* global_MBR = GlobalMbr_parallel(&geometries, rank);

    // Creating grid over MBR in root only
    vector<vector<GEOSGeometry*>> global_MBR_grid = createGrid(ctx, global_MBR, gridSize);

    // Declare a map to store unique ID, geometry, and its sketch
    unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap;

    // Start construction timing for polygon MBR, grids, sketch
    double startTimeSket = MPI_Wtime();
    // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
    for (const auto& pair : centeredGeosWithID) {
        const string& id = pair.first;               // Unique ID of the geometry
        GEOSGeometry* geometry = pair.second;         // Geometry pointer

        vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
        if (sketchResult.empty()) {
            cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
        } else {
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
        std::cout << "total time for sketch creation input after trimming: "
                  << maxTotalTimeSket << " seconds" << std::endl;
    }

    // Creating inverted index
    unordered_map<int, vector<string>> invertedIndex;
    ulong totalVectors = 0;

    for (const auto& pair : geometrySketchMap) {
        const string& id = pair.first;
        const vector<double>& sketch = pair.second.second;
        for (size_t i = 0; i < sketch.size(); ++i) {
            if (sketch[i] > 0.0) {
                invertedIndex[i].push_back(id);
            }
        }
        totalVectors++;
    }

    // Store the length of invertedIndex for each grid cell
    unordered_map<int, int> invertedIndexLength;
    for (int i = 0; i < gridSize * gridSize; ++i) {
        if (invertedIndex.find(i) == invertedIndex.end()) {
            invertedIndexLength[i] = 0;
        } else {
            invertedIndexLength[i] = invertedIndex[i].size();
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Save the inverted index lengths to a file
    string outputFilename = "inverted_index_lengths" + to_string(rank) + ".csv";
    ofstream outFile(
        "/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/inverted_index_lengths/"
        + outputFilename,
        ios::out
    );
    if (!outFile.is_open()) {
        cerr << "Failed to open file for writing inverted index lengths." << endl;
        MPI_Finalize();
    }
    outFile << "index,length\n";
    for (size_t i = 0; i < invertedIndexLength.size(); ++i) {
        outFile << i << "," << invertedIndexLength[i] << endl;
    }
    outFile.close();

    // To store the length of the seen_ids set
    vector<unsigned int> seen_id_Lengths;

    // Start the timing for the query processing
    double startTimeQuery = MPI_Wtime();

    // Now, let's process each query:
    for (const auto& queryPair : centeredQueryGeosWithID) {
        const string& queryID = queryPair.first;    // Query geometry ID
        GEOSGeometry* queryGeometry = queryPair.second;  // Query geometry pointer

        // Generate sketch for the query geometry (trimmed if needed)
        vector<double> querySketch = sketch(ctx, &global_MBR_grid, queryGeometry);
        double querySketchSum = sketchSum(querySketch);
        double traversalSum = querySketchSum;

        // Result vector for storing the top-k neighbors (Jaccard distances and geometry IDs)
        vector<pair<double, string>> queryResults;

        // Vector to store non-zero grid cells (flattened indexes)
        vector<int> nonZeroCellIndexes;

        // Iterate over each non-zero grid cell in the query sketch
        for (size_t i = 0; i < querySketch.size(); ++i) {
            if (querySketch[i] > 0.0) {
                nonZeroCellIndexes.push_back(i);
            }
        }

        // Sort the non-zero cell indexes in ascending order by inverted list size
        sort(nonZeroCellIndexes.begin(), nonZeroCellIndexes.end(),
             [&invertedIndex](int a, int b) {
                 return invertedIndex[a].size() < invertedIndex[b].size();
             });

        // Define the seen_ids set to avoid duplicate processing
        unordered_set<string> seen_ids;

        for (int flattenedIndex : nonZeroCellIndexes) {
            if (invertedIndex.find(flattenedIndex) != invertedIndex.end()) {
                const vector<string>& relevantGeometries = invertedIndex[flattenedIndex];

                // For each geometry in this grid cell, calculate the Jaccard distance
                for (const auto& geoID : relevantGeometries) {
                    if (seen_ids.find(geoID) != seen_ids.end()) {
                        continue;
                    }
                    const auto& geoPair = geometrySketchMap[geoID];
                    GEOSGeometry* geometry = geoPair.first;
                    double jaccardDist = jaccardDistance(ctx, queryGeometry, geometry);
                    queryResults.push_back(make_pair(jaccardDist, geoID));
                    seen_ids.insert(geoID);
                }
            }

            traversalSum -= querySketch[flattenedIndex];
            if (traversalSum <= 0.7 * querySketchSum) {
                break;
            }
        }

        // Sort the results by Jaccard distance (ascending order)
        sort(queryResults.begin(), queryResults.end(),
             [](const pair<double, string>& a, const pair<double, string>& b) {
                 return a.first < b.first;
             });

        // Store the size of the seen_ids set in the seen_id_Lengths vector
        seen_id_Lengths.push_back(seen_ids.size());

        // Keep only the top-k nearest neighbors
        if (queryResults.size() > numNeighbor) {
            queryResults.resize(numNeighbor);
        }

        // Store the query ID and the list of top-k results in PrefixFilteroutput
        PrefixFilteroutput.push_back(make_pair(queryID, queryResults));
    }

    // End the timing for the query processing
    double endTimeQuery = MPI_Wtime();
    double totalTimeQuery = endTimeQuery - startTimeQuery;
    double maxTotalTimeQuery;

    // Use MPI_Reduce to find the maximum query processing time across all processes
    MPI_Reduce(&totalTimeQuery, &maxTotalTimeQuery, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Maximum Query Processing in Spatial Prefix Filter: "
             << maxTotalTimeQuery << " seconds" << endl;
    }

    // Save the seen_id_Lengths to a file
    string outputSeenIDLengthsFilename = "seen_id_lengths" + to_string(rank) + ".csv";
    ofstream outSeenIDLengthsFile(
        "/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/seen_id_lengths/"
        + outputSeenIDLengthsFilename,
        ios::out
    );
    if (!outSeenIDLengthsFile.is_open()) {
        cerr << "Failed to open file for writing seen ID lengths." << endl;
        MPI_Finalize();
    }
    for (size_t i = 0; i < seen_id_Lengths.size(); ++i) {
        outSeenIDLengthsFile << seen_id_Lengths[i];
        if (i != seen_id_Lengths.size() - 1) {
            outSeenIDLengthsFile << ",";
        }
    }
    outSeenIDLengthsFile << "\n";
    outSeenIDLengthsFile.close();
}

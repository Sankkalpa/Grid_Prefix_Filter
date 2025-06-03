
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


int main(int argc, char** argv) {
    int rank, size;
    int numSplit;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Comm comm = MPI_COMM_WORLD;
    char hostname[256];
    gethostname(hostname, 256);

    if (argc != 4) {
        if (rank == 0) {
            cerr << "Incorrect number of arguments!\n";
            cerr << "Usage: mpirun -n <num_processes> ./parallel_main <dataFile> <inputFile>\n";
        }
        MPI_Finalize();
        return 1;
    }

    string datafile  = string(argv[1]);
    string inputfile = string(argv[2]);
    unsigned int numNeighbor = atoi(argv[3]);

    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    GEOSWKTWriter* wktWriter = GEOSWKTWriter_create_r(ctx);

    cout << "error1" << endl;
    vector<string> localFiles = splitFileMPI(datafile, comm);
    cout << "error2" << endl;

    unordered_map<string, GEOSGeometry*>* geos = read_wkt_with_id(localFiles[0], ctx);
    cout << "error3" << endl;

    // Print ID + geometries
    if (geos != nullptr && !geos->empty()) {
        for (const auto& pair : *geos) {
            const string& id       = pair.first;
            GEOSGeometry* geometry = pair.second;

            cout << "ID: " << id << " - input Geometry: ";
            printCoordinates(ctx, geometry);
            cout << endl;
        }
    } else {
        cerr << "No geometries to print or failed to read geometries." << endl;
    }

    if (geos->empty()) {
        cerr << "Error: No geometries to process in process " << rank << endl;
        MPI_Finalize();
        GEOS_finish_r(ctx);
        return 1;
    }

    // Read query geometries and their unique ID
    unordered_map<string, GEOSGeometry*>* querygeos = read_wkt_with_id(inputfile, ctx);

    // Vector to hold centered query geometries
    unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID;
    for (const auto& pair : *querygeos) {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            centeredQueryGeosWithID[pair.first] = centeredGeo;
        } else {
            cerr << "Error centering geometry with ID " << pair.first << endl;
        }
    }

    // Perform prefix-filter query
    int gridSize = 50;
    ResultType PrefixFilteroutput;

    if (rank == 0) {
        cout << "Executing PrefixFilterFunction..." << endl;
    }
    PrefixFilterFunction(
        *geos,
        rank,
        centeredQueryGeosWithID,
        numNeighbor,
        gridSize,
        PrefixFilteroutput,
        ctx
    );
    if (rank == 0) {
        cout << "PrefixFilterFunction completed." << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Write prefix-filter results to file
    string outputFilename = "Prefix_results_rank_" + to_string(rank) + ".txt";
    ofstream outFile("path_to_results/Prefix_results_rank_" + to_string(rank) + ".txt", ios::out);
    if (!outFile.is_open()) {
        cerr << "Failed to open file for writing results." << endl;
        MPI_Finalize();
        return 1;
    }
    for (const auto& queryResult : PrefixFilteroutput) {
        const string& queryID   = queryResult.first;
        const auto& neighbors   = queryResult.second;

        outFile << "Query ID: " << queryID << endl;
        for (const auto& neighbor : neighbors) {
            outFile << "Neighbor ID: " << neighbor.second
                    << ", Distance: " << neighbor.first << endl;
        }
        outFile << endl;
    }
    outFile.close();

    MPI_Barrier(MPI_COMM_WORLD);

    // Aggregate prefix-filter results on rank 0
    unordered_map<string, QueryResult> globalResults;
    if (rank == 0) {
        string baseDirectory  = "path_to_results/results/";
        string baseFilename   = "Prefix_results_rank_";
        string fileExtension  = ".txt";
        unordered_map<string, QueryResult> allResult;

        for (int i = 0; i < size; ++i) {
            string filename = baseDirectory + baseFilename + to_string(i) + fileExtension;
            ifstream inFile(filename);

            if (inFile.is_open()) {
                allResult.clear();
                ParseFromFile(ctx, filename, allResult);

                for (const auto& [queryID, neighbors] : allResult) {
                    vector<NeighborPair> combinedNeighbors = globalResults[queryID];
                    combinedNeighbors.insert(combinedNeighbors.end(),
                                             neighbors.begin(),
                                             neighbors.end());
                    sort(combinedNeighbors.begin(), combinedNeighbors.end(),
                         [](const NeighborPair& a, const NeighborPair& b) {
                             return a.first < b.first;
                         });
                    if (combinedNeighbors.size() > numNeighbor) {
                        combinedNeighbors.resize(numNeighbor);
                    }
                    globalResults[queryID] = combinedNeighbors;
                }
                inFile.close();
            } else {
                cerr << "Failed to open file: " << filename << endl;
            }
        }
    }

    // Brute-force query
    vector<pair<string, vector<pair<double, string>>>> queryBFResults;
    BFquery(rank, centeredQueryGeosWithID, *geos, numNeighbor, queryBFResults);

    MPI_Barrier(MPI_COMM_WORLD);

    // Write brute-force results to file
    string outputBFFilename = "BFquery_results_rank_" + to_string(rank) + ".txt";
    ofstream outBFFile(
        "path_to_file/src/data/results/" + outputBFFilename,
        ios::out
    );
    if (!outBFFile.is_open()) {
        cerr << "Failed to open file for writing results." << endl;
        MPI_Finalize();
        return 1;
    }
    for (const auto& resultPair : queryBFResults) {
        const string& queryID         = resultPair.first;
        const auto& neighborPairs     = resultPair.second;

        outBFFile << "Query ID: " << queryID << endl;
        for (const auto& neighborPair : neighborPairs) {
            double distance          = neighborPair.first;
            const string& neighborID = neighborPair.second;
            outBFFile << "Neighbor ID: " << neighborID
                      << ", Distance: " << distance << endl;
        }
        outBFFile << endl;
    }
    outBFFile.close();

    MPI_Barrier(MPI_COMM_WORLD);

    // Aggregate brute-force results on rank 0
    unordered_map<string, QueryResult> globalBFResults;
    if (rank == 0) {
        string baseDirectory  = "path_to_file/src/data/results/";
        string baseFilename   = "BFquery_results_rank_";
        string fileExtension  = ".txt";
        unordered_map<string, QueryResult> allBFResult;

        for (int i = 0; i < size; ++i) {
            string filename = baseDirectory + baseFilename + to_string(i) + fileExtension;
            ifstream inFile(filename);

            if (inFile.is_open()) {
                allBFResult.clear();
                ParseFromFile(ctx, filename, allBFResult);

                for (const auto& [queryID, neighbors] : allBFResult) {
                    vector<NeighborPair> combinedNeighbors = globalBFResults[queryID];
                    combinedNeighbors.insert(combinedNeighbors.end(),
                                             neighbors.begin(),
                                             neighbors.end());
                    sort(combinedNeighbors.begin(), combinedNeighbors.end(),
                         [](const NeighborPair& a, const NeighborPair& b) {
                             return a.first < b.first;
                         });
                    if (combinedNeighbors.size() > numNeighbor) {
                        combinedNeighbors.resize(numNeighbor);
                    }
                    globalBFResults[queryID] = combinedNeighbors;
                }
                inFile.close();
            } else {
                cerr << "Failed to open file: " << filename << endl;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Compare performance on rank 0
    if (rank == 0) {
        FILE* output = fopen(
            "path_to_file/src/data/comparisonResults.csv",
            "w"
        );
        if (output == NULL) {
            cerr << "Failed to open file for writing.\n";
        } else {
            fprintf(output,
                    "query_id,avg_recall,mse,false_positive_rate,"
                    "total_lsh_neighbors,total_bf_neighbors,total_matched_neighbors\n");

            double totalOverallRecall        = 0.0;
            double totalOverallBFRecall      = 0.0;
            double totalMSE                  = 0.0;
            int totalCount                   = 0;
            int totalFalsePositives          = 0;
            int totalBFNeighborCount         = 0;
            int totalLSHNeighborCount        = 0;
            int totalMatchedNeighborsCount   = 0;

            for (const auto& [queryID, bfNeighbors] : globalBFResults) {
                auto lshResultIt = globalResults.find(queryID);
                if (lshResultIt == globalResults.end()) {
                    continue;
                }
                auto& lshNeighbors = lshResultIt->second;
                int correctNeighbors = 0;
                double mse        = 0.0;
                int mseidx        = 0;

                totalBFNeighborCount  += bfNeighbors.size();
                totalLSHNeighborCount += lshNeighbors.size();

                for (const auto& [distance, lshNeighborID] : lshNeighbors) {
                    bool isMatched = false;
                    for (const auto& [bfDistance, bfNeighborID] : bfNeighbors) {
                        if (lshNeighborID == bfNeighborID) {
                            correctNeighbors++;
                            isMatched = true;
                            mse += pow(bfDistance - distance, 2);
                            mseidx++;
                            break;
                        }
                    }
                    if (!isMatched) {
                        totalFalsePositives++;
                    }
                }

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

            double averageRecall       = (totalCount > 0) ? totalOverallRecall / totalCount : 0;
            double averageBFRecall     = (totalCount > 0) ? totalOverallBFRecall / totalCount : 0;
            double averageMSE          = (totalCount > 0) ? totalMSE / totalCount : 0;
            double falsePositiveRate   = (totalLSHNeighborCount > 0)
                                            ? static_cast<double>(totalFalsePositives) / totalLSHNeighborCount
                                            : 0;

            cout << "Average Recall: " << averageRecall << endl;
            cout << "Average Recall by BF neighbors: " << averageBFRecall << endl;
            cout << "MSE: " << averageMSE << endl;
            cout << "Total Matched Neighbors: " << totalMatchedNeighborsCount << endl;
            cout << "Total Neighbors (including duplicates): " << totalLSHNeighborCount << endl;
            cout << "Total BF Neighbors (including duplicates): " << totalBFNeighborCount << endl;
            cout << "Total False Positives rate: " << falsePositiveRate << endl;

            fclose(output);
        }
    }

    MPI_Finalize();
    GEOS_finish_r(ctx);
    return 0;
}



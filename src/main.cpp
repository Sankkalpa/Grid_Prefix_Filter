
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



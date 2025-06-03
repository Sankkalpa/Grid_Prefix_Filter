#include "query.h"
#include <mpi.h>

using namespace std;



void BFquery(
    int rank,
    const unordered_map<string, GEOSGeometry*>& queryGeos,
    const unordered_map<string, GEOSGeometry*>& geos,
    unsigned int numNeighbors,
    vector<pair<string, vector<pair<double, string>>>>& queryResults
) {
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    unordered_map<string, GEOSGeometry*> centeredGeosWithID;
    // Center each geometry and store it with its identifier
    for (const auto& pair : geos) {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            centeredGeosWithID[pair.first] = centeredGeo;
        } else {
            std::cerr << "Error centering geometry with ID "
                      << pair.first << std::endl;
        }
    }

    // Start timing BF
    double startTimeBF = MPI_Wtime();

    for (const auto& queryPair : queryGeos) {
        const string& queryID = queryPair.first;
        GEOSGeometry* queryGeo = queryPair.second;
        vector<pair<double, string>> nn;  // Nearest neighbors

        for (const auto& candidatePair : centeredGeosWithID) {
            const string& candidateID = candidatePair.first;
            GEOSGeometry* candidateGeo = candidatePair.second;
            double distance = jaccardDistance(ctx, queryGeo, candidateGeo);
            nn.emplace_back(distance, candidateID);
        }

        // Sort the nearest neighbors by distance
        sort(nn.begin(), nn.end());
        // Keep only the top numNeighbors neighbors
        nn.resize(min(nn.size(), static_cast<size_t>(numNeighbors)));

        queryResults.push_back(make_pair(queryID, nn));
    }

    double endTimeBF = MPI_Wtime();
    double totalTimeBF = endTimeBF - startTimeBF;
    double maxTimeBF;

    // Use MPI_Reduce to find the maximum BF query time across all processes
    MPI_Reduce(&totalTimeBF, &maxTimeBF, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Print the maximum BF query time in the master process (rank 0)
    if (rank == 0) {
        std::cout << "Maximum BF query time from all processes: "
                  << maxTimeBF << " seconds" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    GEOS_finish_r(ctx);
}

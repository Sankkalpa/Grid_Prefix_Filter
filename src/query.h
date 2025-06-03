// Query types

#ifndef QUERY_H
#define QUERY_H

#include "geoutil.h"
#include "thread.h"
#include "util.h"
#include <set>
#include <mpi.h>
#include <unordered_set>

// LSH queries

using ResultType = std::vector<std::pair<std::string, std::vector<std::pair<double, std::string>>>>;

double sketchSum(const std::vector<double>& sketch);

void BFquery(
             int rank,
             const unordered_map<string, GEOSGeometry*>& queryGeos,
             const unordered_map<string, GEOSGeometry*>& geos,
             unsigned int numNeighbors,
             vector<pair<string, vector<pair<double, string>>>>& queryResults
             ) ;

void PrefixFilterFunction(
             const unordered_map<string, GEOSGeometry*>& geos, // Changed to unordered_map
             int rank,
             const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
             unsigned int numNeighbor,
             int gridSize,
             ResultType& PrefixFilteroutput,
             GEOSContextHandle_t ctx
            ) ;



#endif


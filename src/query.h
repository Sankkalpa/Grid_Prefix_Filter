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
             
void queryThreadFunct(
             unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
             unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap,
             vector<HashMap<std::string>*> hashMaps,
             int numNeighbors, 
             int hashLength, 
             vector<unsigned int>& seeds, 
             int rank,
             ResultType& LSHoutput,
             vector<vector<GEOSGeometry*>>& global_MBR_grid,
             const vector<bool>& cellsToKeep
             ) ;
                     
void LSHHashFunction(
             const unordered_map<string, GEOSGeometry*>& geos, // Changed to unordered_map
             int rank,
             const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
             unsigned int numNeighbor,
             int gridSize,
             int hashLength ,
            // std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>>& LSHoutput,
             ResultType& LSHoutput,
             GEOSContextHandle_t ctx
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


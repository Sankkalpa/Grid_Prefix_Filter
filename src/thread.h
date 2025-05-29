// Definition of thread body functions as well as input and output structures.
// Some of the structures have very little inside them, but are still wrapped
// in structures to maintain consistency.

#ifndef THREAD_H
#define THREAD_H

#define GEOS_USE_ONLY_R_API
#include "geos_c.h"
#include "geoutil.h"
#include "parse_geodata.h"
#include "util.h"
#include <algorithm>
#include <barrier>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>

// ----------------------------------------------------------------------------
// Input Threads : Read data from file and preprocess it by centering the
//                 polygons and tracking a minimum bounding rectangle.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TYPE inputIn :
//     Bundled input for the input threads to ease with managment of input to
//     these threads.
// MEMBERS :
//     string filename : Name of the file to read data from.
// ----------------------------------------------------------------------------
typedef struct inputIn
{
    string filename; // Name of file to read from
} inputIn;

// ----------------------------------------------------------------------------
// TYPE inputOut :
//     Bundled output for the input threads.
// MEMBERS :
//     vector<GEOSGeometry *> *geos : All of the processed polygons.
// ----------------------------------------------------------------------------
typedef struct inputOut
{
    vector<GEOSGeometry *> *geos;
} inputOut;

// ----------------------------------------------------------------------------
// FUNCTION inputThreadFunc :
//     Body function for the input threads. Reads data from file, center it
//     and track a local minimum bounding rectangle.
// PARAMETERS :
//     inputIn *in   : Collective input for the thread.
//     inputOut *out : Collective output for the thread.
// RETURNS : Technically void, all generated data stored in OUT.
// ----------------------------------------------------------------------------
void inputThreadFunc(inputIn *in, inputOut *out);

// ----------------------------------------------------------------------------
// Construction Threads : Create sketches of the polygons, then generate the
//                        LSH hashes and insert them into the hash maps.
// ----------------------------------------------------------------------------

/* NOTE: THE ORDER OF THE NEXT THREE DECLARTIONS IS VERY SPECIFIC ON PURPOSE */

// ----------------------------------------------------------------------------
// TYPE constructionOut :
//     Bundled output for the construction threads to ease with managment of
//     parameters within the main thread.
// MEMBERS :
//     vector<int> *cellUsage : Vector to track the usage of each cell of the
//                              grid. Used by the GRIDFILTERER to determine
//                              which cells to drop.
// ----------------------------------------------------------------------------
typedef struct constructionOut
{
    vector<int> *cellUsage;
} constructionOut;

// ----------------------------------------------------------------------------
// CLASS gridFilterer :
//      This class is a bit of an oddity. It is used as a template parameter to
//      the barrier and acts as the completion step of the barrier. Basically,
//      all of the construction threads generate sketches for their polygons
//      but at some point we have to drop grid cells which do not meet usage
//      requirements, which is what this class does. Since I cannot pass any
//      parameters directly to the filter function, the constructor stores
//      references to all needed data even though those pointers point to empty
//      containers until the construction threads generate their output.
// MEMBER VARIABLES :
//     GEOSContextHandle_t ctx             : Thread dependent context handle.
//     vector<vector<GEOSGeometry *> *grid : Grid to filter.
//     vector<constructionOut> *conOut     : Cell usage data from Construction
//                                           threads.
//     list<int> *remove                   : List of all indices which where
//                                           removed.
//     int threshold                       : Number of cell usages required to
//                                           not be deleted.
// METHODS :
//     gridFilterer() : Constructor, just stores parameters to member variables
//     operator()()   : Overloads the () operator, which is what the barrier
//                      object calls upon completion of a phase.
// ----------------------------------------------------------------------------
class gridFilterer
{
  public:
    gridFilterer(GEOSContextHandle_t c, vector<vector<GEOSGeometry *>> *g,
                 vector<constructionOut> *co, int th, list<int> *r)
        : ctx(c), grid(g), conOut(co), remove(r), threshold(th){};
    gridFilterer operator()();

  private:
    GEOSContextHandle_t ctx;
    vector<vector<GEOSGeometry *>> *grid;
    vector<constructionOut> *conOut;
    list<int> *remove;
    int threshold;
};

// ----------------------------------------------------------------------------
// TYPE constructionIn :
//     Bundled input to the construction threads, used to ease parameter
//     management in the main thread.
// MEMBERS :
//     vector<GEOSGeometry *> *geos    : The set of polygons to process.
//     vector<vector<GEOSGeometry *>>
//                               *grid : The grid used for sketching.
//     vector<HashMap<
//                GEOSGeometry *> *>
//                               *maps : The set of hash maps to insert each
//                                       polygon into.
//     int hashLenght                  : Number of elements in each LSH hash.
//     vector<vector<unsigned int>>
//                              *seeds : Random seeds used to generate each
//                                       minhash within the final LSH hash.
//     barrier<gridFilterer> *bar      : Barrier to wait at for cell usage
//                                       processing by the GRIDFILTERER.
//     list<int> *removedCells         : Output of BAR. This list is empty
//                                       until the return from waiting on the
//                                       barrier.
// ----------------------------------------------------------------------------
typedef struct constructionIn
{
    vector<GEOSGeometry *> *geos;
    vector<vector<GEOSGeometry *>> *grid;
    vector<HashMap<GEOSGeometry *> *> *maps;
    int hashLength;
    vector<vector<unsigned int>> *seeds;
    barrier<gridFilterer> *bar;
    list<int> *removedCells;
} constructionIn;

// ----------------------------------------------------------------------------
// FUNCTION constructionThreadFunc :
//     Body function for the construction thread. Generate, then refine the
//     sketches for a set of polygons before inputing them to the hash maps.
//
//     Each polygon has a different hash generated for each hash map and is
//     inserted into each hash map in order to increase the chances of finding
//     true nearest neighbors.
// PARAMETERS :
//     constructionIn *in   : Collective input to a construction thread.
//     constructionOut *out : Collective output from a construction thread.
// RETURNS : Technically void, all output is placed in OUT.
// ----------------------------------------------------------------------------
void constructionThreadFunc(constructionIn *in, constructionOut *out);

// ----------------------------------------------------------------------------
// Query Threads : Read data from a query file, generate a sketch and then set
//                 of LSH threads for each polygon and linearlly scan the
//                 union of the buckets to find approximate nearest neighbors.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TYPE queryIn :
//     Bundled input to a query thread. Used to ease parameter management in
//     the main thread.
// MEMBERS :
//     string filename                 : The relative path to the query file.
//     unsigned int numNeighbor        : The number of approximate nearest
//                                       neighbors to find.
//     vector<HashMap<
//                GEOSGeometry *> *>
//                               *maps : The set of hash maps to search.
//     int hashLenght                  : Number of elements in each hash.
//     vector<vector<GEOSGeometry *>>
//                               *grid : Grid used to generate sketches of the
//                                       read input polygons.
//     vector<vector<unsigned int>>
//                              *seeds : Random seeds used to generate proper
//                                       LSH hashes in a consitant manner as
//                                       the construction threads.
// ----------------------------------------------------------------------------
typedef struct queryIn
{
    vector<GEOSGeometry *> *queryGeos;
    unsigned int numNeighbor;
    vector<HashMap<GEOSGeometry *> *> *maps;
    int hashLength;
    vector<vector<GEOSGeometry *>> *grid;
    vector<vector<unsigned int>> *seeds;
} queryIn;

// ----------------------------------------------------------------------------
// TYPE queryOut :
//     Bundled output of the query threads. Used to ease parameter management
//     in the main thread.
// MEMBERS :
//     list<pair<GEOSGeometry *,
//         vector<pair<double, GEOSGeometry *>>>>
//                                      *results : The results of the query.
//                                                 Data is structured as a list
//                                                 of pairs where each pair is
//                                                 the result of querying for
//                                                 one of the input polygons.
//                                                 The outer pair maps that
//                                                 GEOSGeometry* to a sorted
//                                                 vector of approximate
//                                                 nearest neighbors, which
//                                                 itself maps the Jaccard
//                                                 Distance bewteen the query
//                                                 shape and the neighbor
//                                                 itself. For convience, this
//                                                 type as been aliased to
//                                                 "QueryResultSet"
// ----------------------------------------------------------------------------
typedef struct queryOut
{
    QueryResultSet *results;
} queryOut;

// ----------------------------------------------------------------------------
// FUNCTION queryThreadFunc :
//     Body function for a query thread. Conducts the entire query process for
//     a set of input polygons. Generate sketchs, generate LSH hashes and scan
//     the union of the buckets from all of the hash maps.
// PARAMETERS :
//     queryIn *in   : Collective input for the query thread.
//     queryOut *out : Collective output for the query thread.
// RETURNS : Technically void, all output placed in the OUT structure.
// ----------------------------------------------------------------------------
void queryThreadFunc(queryIn *in, queryOut *out);
#endif

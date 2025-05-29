// Header file for parsing data from UCR Star datasets
// All of this functions return an array of GEOSGeometry pointers.
// It is the callers responsibilty to free each of these pointers
// with GEOSGeom_destory once they are no longer needed.

#ifndef GEO_PARSER
#define GEO_PARSER

#define GEOS_USE_ONLY_R_API
//#include "geos_c.h"
#include "util.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
#include <unordered_map>


#define INIT_INPUT_SIZE 10

using namespace std;

// ----------------------------------------------------------------------------
// FUNCTION read_csv :
//     Read a CSV file which has wkt as the first column. Uses an outdated
//     return type but it's used in the rest of the project.
// PARAMETERS :
//     const char *filename    : Relative path to the input file.
//     unsigned int *numGeo    : Place to return the number of created
//                               polygons.
//     GEOSContextHandle_t ctx : Thread dependent context handle.
// RETURNS : Double pointer to the newly created polygons. Should be updated to
//           return a vector<GEOSGeometry*>.
// ----------------------------------------------------------------------------
GEOSGeometry **read_csv(const char *filename, unsigned int *numGeo,
                        GEOSContextHandle_t ctx);

// ----------------------------------------------------------------------------
// FUNCTION read_wkt :
//     Read a WKT file which has WKT embedded in a bit of JSON. If this
//     function causes an error, check if the input file's first line is
//     missing a space at the very being. This function makes assumptions about
//     where the WKT starts within the JSON based on the number of character
//     into the line.
// PARAMETERS :
//     const string filename   : Relative path to the input file.
//     GEOSContextHandle_t ctx : Thread dependent context handle.
// RETURNS : Vector of the internal GEOS representation of the shapes in
//           FILENAME, filtered to be only the polygons and it's unique id's.
// ----------------------------------------------------------------------------
unordered_map<string, GEOSGeometry*> *read_wkt_with_id(const string &filename, GEOSContextHandle_t ctx);



                                          
#endif

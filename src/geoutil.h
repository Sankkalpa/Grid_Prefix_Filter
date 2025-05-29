// geometric utilties using GEOS
//
// These are small helper functions which interact directly with GEOS shapes but
// are otherwise utilities.

#ifndef GEOUTIL_H
#define GEOUTIL_H

#define GEOS_USE_ONLY_R_API
#include "geos_c.h"
#include "util.h"
#include <unordered_map>
#include <string>




// ----------------------------------------------------------------------------
// FUNCTION centerGeometry :
//     This function translates a geometry so that it's center is at (0,0).
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle
//     GEOSGeometry *g         : The geometry to be translated
// RETURNS : Pointer to translated geometry
// ----------------------------------------------------------------------------
GEOSGeometry *centerGeometry(GEOSContextHandle_t ctx, GEOSGeometry *g);

// ----------------------------------------------------------------------------
// FUNCTION createGrid :
//     If the BASE is a rectangle, tile it with smaller rectangles so that it
//     creates a GRIDSIZE x GRIDSIZE grid over the same area as base.
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle
//     GEOSGeometry *base      : Rectangle whic is to be divided
//     int gridSize            : Number of cells which should line each side of
//                               the generated grid.
// RETURNS : Generated grid as a 2D array (vector of vectors)
// ----------------------------------------------------------------------------
vector<vector<GEOSGeometry *>> createGrid(GEOSContextHandle_t ctx,
                                          GEOSGeometry *base, int gridSize);

// ----------------------------------------------------------------------------
// FUNCTION fillRatio :
//     Determine the proportion of BASE's area which is also within OVERLAY.
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle
//     GEOSGeometry *base      : The polygon we are interesting in "filling"
//     GEOSGeometry *overlay   : The polygon which does the filling.
// RETURNS : Proportion of BASE's area also in OVERLAY.
// ----------------------------------------------------------------------------
double fillRatio(GEOSContextHandle_t ctx, GEOSGeometry *base,
                 GEOSGeometry *overlay);

// ----------------------------------------------------------------------------
// FUNCTION getCenter :
//     Find the centroid of a polgyon.
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle
//     GEOSGeometry *g         : Geometry to find the center of.
//     double *x               : Where to store the x coordinate.
//     double *y               : Where to store the y coordinate.
// RETURNS : Technically void, but output is stored in double* passed in.
// ----------------------------------------------------------------------------
void getCenter(GEOSContextHandle_t ctx, GEOSGeometry *g, double *x, double *y);

// ----------------------------------------------------------------------------
// FUNCTION jaccardDistance :
//     Find the Jaccard Distance between to shapes. Defined as the one minus
//     the ratio of intersecting area over the area of the union of the shapes.
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle.
//     GEOSGeometry *g1        : First of the two shapes.
//     GEOSGeometry *g2        : Second of the two shapes.
// RETURNS : Double representing the Jaccard distance between the shapes.
// ----------------------------------------------------------------------------
double jaccardDistance(GEOSContextHandle_t ctx, GEOSGeometry *g1,
                       GEOSGeometry *g2);
                       
                       
                       
// ----------------------------------------------------------------------------
// FUNCTION jaccardDistance :
//     Find the Jaccard Distance between to sketch. Defined as the one minus
//     the ratio of minimum area over the maximum area the shapes sketch.
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle.
//     std::vector<double>& sketch1: First of the two sketch.
//     std::vector<double>& sketch2: Second of the two sketch.
// RETURNS : Double representing the Jaccard distance between the sketch.
// ----------------------------------------------------------------------------                       
double SketchJaccardDistance(const std::vector<double>& sketch1, const std::vector<double>& sketch2) ;                      

// ----------------------------------------------------------------------------
// FUNCTION LSHHash :
//     Compute a locality sensitive hash using the method propoused by
//     A. Shrivastava in "Simple and Efficient Weighted Minwise Hashing".
//
//     This algorithm operates on vectors and divides the numberline into
//     regions for each element of the vector based on the maximum value
//     observed within the dataset for that element. The green region starts at
//     bottom of each element's region and extends upwards for the value that
//     this vector has in that element. The rest of the area is the red region
//
//     It then generates random numbers over this total range and generates a
//     minhash of the number of attempts it took to generate a number in the
//     green region. Each minhash becomes one element in the final hash.
// PARAMETERS :
//     vector<double> *sketch      : The sketch of the polygon (see below).
//     int hashLength              : The hash is repsented as a vector, so we
//                                   can generate a hash of arbitary length.
//                                   This determines the length of the hash.
//     vector<unsigned int> *seeds : Each element of the hash needs a random
//                                   seed to start the random number generator
//                                   so that the process is consistant over the
//                                   whole dataset.
// RETURNS : A locality senstive hash of the input SKETCH.
// ----------------------------------------------------------------------------
vector<int> LSHHash(vector<double> *sketch, int hashLength,
                    vector<unsigned int> *seeds);
// ----------------------------------------------------------------------------
// FUNCTION LSHHash :
//Same as above LSHHash function but only takes one seed for entire hashing at a time
// ----------------------------------------------------------------------------                    
                    
vector<int> LSHHashs(vector<double> *sketch, int hashLength, unsigned int seed);                    

// ----------------------------------------------------------------------------
// FUNCTION minimumBoundingRectangle :
//     Find the minimum bounding rectangle over a set of input geometries.
// PARAMETERS :
//     GEOSContextHandle_t ctx    : Thread dependent context handle.
//     vector<GEOSGeometry *> geo : Set of input geometries.
// RETURNS : GEOSGeometry* repersenting the minimum bounding rectangle of GEO
// ----------------------------------------------------------------------------

GEOSGeometry *minimumBoundingRectangle(GEOSContextHandle_t ctx, GEOSGeometry *geometry);
// ----------------------------------------------------------------------------
// FUNCTION sketch:
//     Generate a sketch of a polygon. To do this, we require a grid over the
//     entire area of interest. Consider one cell within that grid. It is
//     filled to some level (possibly zero) by the polygon we are created the
//     sketch of. Each element of the output vector corresponds to one such
//     cell of the input grid and has a value equal to the proportion of the
//     grid's area which is filled by the input polygon.
//
//     NOTE: Some grid elements may be NULL, in which case they are skipped.
//     This is an appoarch to reduce processing time with limited loss of
//     accuracy.
// PARAMETERS :
//     GEOSContextHandle_t ctx            : Thread dependent context handle.
//     vector<vector<GEOSGeometry *> grid : Grid over global area of interest.
//     GEOSgeometry *g                    : Polygon to generate a sketch over.
// RETURNS : A vectoized representation of the input polygon.
// ----------------------------------------------------------------------------
vector<double> sketch(GEOSContextHandle_t ctx,
                      vector<vector<GEOSGeometry *>> *grid, GEOSGeometry *g);
                      
// ----------------------------------------------------------------------------
// FUNCTION printCoordinates:
//prints the coordinates of the polygons
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle.
//     GEOSGeometry *geometry  : Polygon geometery to be printed

// RETURNS : coordinates of the input polygon.
// ---------------------------------------------------------------------------- 
void printCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry);
                     
 
unordered_map<string, GEOSGeometry*> copyGeometriesMap(const unordered_map<string, GEOSGeometry*>& originalMap, GEOSContextHandle_t ctx);

//to find bounding box for a polygon
std::pair<std::pair<double, double>, std::pair<double, double>> findBoundingBox(const std::vector<std::pair<double, double>>& points);

//to calculate local mbr's for each processes
GEOSGeometry *calculateAggregateLocalMBR(GEOSContextHandle_t ctx, const vector<GEOSGeometry *> &localMBRs);

// Function to calculate the envelope of a geometry
GEOSGeometry *calculateEnvelope(GEOSContextHandle_t ctx, GEOSGeometry *geometry);

std::vector<std::pair<double, double>> extractCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry);

//to create rectangle over boundingbox coordinates
GEOSGeometry *createPolygonFromBoundingBox(GEOSContextHandle_t ctx, const std::pair<std::pair<double, double>, std::pair<double, double>> &boundingBox);
 
//to calculate global MBR at root 
GEOSGeometry  *GlobalMbr_parallel(vector<GEOSGeometry *> *geos, int rank);            
#endif

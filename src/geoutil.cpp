#include "geoutil.h"
#include <mpi.h>

vector<GEOSGeometry *> rootAggregateMBRs;
GEOSGeometry *boundingBoxGeometry = nullptr;
GEOSGeometry *globalMBR = nullptr;


//to gather MBR coordinates at root from each processes
std::vector<std::pair<double, double>> gatherCoordinatePairs(const std::vector<std::pair<double, double>>& coordinates, int rank)
{
    std::vector<std::pair<double, double>> allCoordinates;

    // Calculate the total number of processes
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate space on the root process to store gathered values
    if (rank == 0)
    {
        allCoordinates.resize(size * coordinates.size());
    }

    // Gather values from all processes to the root process
    MPI_Gather(coordinates.data(), 2 * static_cast<int>(coordinates.size()), MPI_DOUBLE,
               allCoordinates.data(), 2 * static_cast<int>(coordinates.size()), MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Return the gathered coordinates to all processes
    return allCoordinates;
}




//To create Global MBR at root and then broadcast it to every other processes
GEOSGeometry  *GlobalMbr_parallel(vector<GEOSGeometry *> *geos, int rank)
{ 
// Create new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    
    GEOSGeometry *aggregateLocalMBR = nullptr;
    vector<GEOSGeometry *> localMBRs;

    // Calculate MBRs for each polygon in the specified range
   for (std::vector<GEOSGeom_t*>::size_type i = 0; i < geos->size(); ++i)
    {
        GEOSGeometry *currentPolygon = geos->at(i);
      
        GEOSGeometry *mbr = minimumBoundingRectangle(ctx, currentPolygon);
        if (!mbr)
        {
            cerr << "Error: Null MBR for polygon " << i << " in process " << rank << endl;
            continue;
        }
        localMBRs.push_back(mbr); // Store the local MBR in the vector
        
        GEOSGeom_destroy_r(ctx, currentPolygon);
       
    }

    // Calculate aggregate local MBR for each process
    GEOSGeometry *LocalMBR_union = calculateAggregateLocalMBR(ctx, localMBRs);
     
    
    aggregateLocalMBR = minimumBoundingRectangle(ctx, LocalMBR_union);
    
    GEOSGeom_destroy_r(ctx, LocalMBR_union);
   
    //convert mbr coordinates to double type 
    std::vector<std::pair<double, double>> mbrCoordinates = extractCoordinates(ctx, aggregateLocalMBR);
  
    //gather all double coordinates to root
    std::vector<std::pair<double, double>> mbrCoordinates_gathered  = gatherCoordinatePairs(mbrCoordinates,  rank);
      
    //extracting only MBR coordinates
     auto boundingBox = findBoundingBox(mbrCoordinates_gathered);
     
    //broadcasting double type MBR coordinates to all process
    MPI_Bcast(&boundingBox, sizeof(std::pair<std::pair<double, double>, std::pair<double, double>>), MPI_BYTE, 0, MPI_COMM_WORLD);
 
    //creating geos polygon from double data type
    boundingBoxGeometry = createPolygonFromBoundingBox(ctx, boundingBox);
     
    globalMBR = minimumBoundingRectangle(ctx, boundingBoxGeometry);
   
    GEOSGeom_destroy_r(ctx, boundingBoxGeometry);

return globalMBR;

}
    

//to copy geometries with it's unique id 
unordered_map<string, GEOSGeometry*> copyGeometriesMap(
const unordered_map<string, GEOSGeometry*>& originalMap, GEOSContextHandle_t ctx
) 
{
    unordered_map<string, GEOSGeometry*> copiedMap;
    for (const auto& pair : originalMap) {
        const string& id = pair.first;
        GEOSGeometry* originalGeometry = pair.second;
        // Clone the geometry using GEOSGeom_clone_r
        GEOSGeometry* copiedGeometry = GEOSGeom_clone_r(ctx, originalGeometry);
        // Insert the id and copied geometry into the new map
        copiedMap.insert({id, copiedGeometry});
    }
    return copiedMap;
}


//to calculate Jaccard distance for comparison
double jaccardDistance(GEOSContextHandle_t ctx, GEOSGeometry *g1,
                       GEOSGeometry *g2)
{
    GEOSGeometry *Intersect = GEOSIntersection_r(ctx, g1, g2);
    double intersectArea;
    GEOSArea_r(ctx, Intersect, &intersectArea);
    GEOSGeom_destroy_r(ctx, Intersect);

    GEOSGeometry *Union = GEOSUnion_r(ctx, g1, g2);
    double unionArea;
    GEOSArea_r(ctx, Union, &unionArea);
    GEOSGeom_destroy_r(ctx, Union);

    return 1.0 - (intersectArea / unionArea);
}

double SketchJaccardDistance(const std::vector<double>& sketch1, const std::vector<double>& sketch2) {
    if (sketch1.size() != sketch2.size()) {
        throw std::invalid_argument("Sketches must be of the same size.");
    }

    double minSum = 0.0; // Sum of minimum values for each cell
    double maxSum = 0.0; // Sum of maximum values for each cell

    for (size_t i = 0; i < sketch1.size(); ++i) {
        minSum += std::min(sketch1[i], sketch2[i]);
        maxSum += std::max(sketch1[i], sketch2[i]);
    }

    // Avoid division by zero in case maxSum is zero
    if (maxSum == 0) {
        return 1.0; // If both sketches are empty, consider them completely dissimilar
    }

    double jaccardSimilarity = minSum / maxSum;

    // Jaccard distance is 1 minus Jaccard similarity
    return 1.0 - jaccardSimilarity;
}

double fillRatio(GEOSContextHandle_t ctx, GEOSGeometry *base,
                 GEOSGeometry *overlay)
{  
 //check for validity
  int isValidb = GEOSisValid_r(ctx, base);

//   cout << "Error 40" << endl;
  
    if (GEOSisValid_r(ctx, base) != 1)
        aprintf(28, "Invalid Grid Cell detected!\n");

// cout << "Error 41" << endl;
     
     int isValido = GEOSisValid_r(ctx, overlay);

    // cout << "Error 42" << endl;
    
    if (GEOSisValid_r(ctx, overlay) != 1) 
        aprintf(32, "Invalid input polygon detected!\n");

    // cout << "Error 43" << endl;
    
    // Only take the intersection if we know that they touch.
    double baseArea;
    GEOSArea_r(ctx, base, &baseArea);

    // cout << "Error 44" << endl;

    GEOSGeometry *intersect = GEOSIntersection_r(ctx, base, overlay);

    // cout << "Error 45" << endl;

    double intersectArea;
    GEOSArea_r(ctx, intersect, &intersectArea);

    // cout << "Error 46" << endl;

    GEOSGeom_destroy_r(ctx, intersect);

    // cout << "Error 47" << endl;

    double result = intersectArea / baseArea;
    if (isnan(result))
        printf("%.4lf / %.4lf = %.4lf\n", intersectArea, baseArea, result);
    return result;
}


void getCenter(GEOSContextHandle_t ctx, GEOSGeometry *g, double *x, double *y)
{
    GEOSGeometry *center = GEOSGetCentroid_r(ctx, g);
    if (GEOSGetNumCoordinates_r(ctx, center) != 1)
    {
        aprintf(39, "Incorrect number of points in centroid\n");
        return;
    }

    GEOSCoordSeq_getXY_r(ctx, GEOSGeom_getCoordSeq_r(ctx, center), 0, x, y);
    GEOSGeom_destroy_r(ctx, center);
}


int transformer(double *x, double *y, void *userdata)
{
    double *centroid = (double *)userdata;
    *x = (*x - centroid[0]) * 100;
    *y = (*y - centroid[1]) * 100;

    return 1;
}


GEOSGeometry *centerGeometry(GEOSContextHandle_t ctx,  GEOSGeometry *g) {
    // Clone the original geometry to ensure it remains unchanged
    GEOSGeometry* gClone = GEOSGeom_clone_r(ctx, g);

    double centroid[2];
    // Use the clone for further operations
    getCenter(ctx, gClone, &(centroid[0]), &(centroid[1]));

    GEOSGeometry *newGeo =
        GEOSGeom_transformXY_r(ctx, gClone, transformer, centroid);
    
    // Destroy the clone after use
    GEOSGeom_destroy_r(ctx, gClone);

    return newGeo; // Return the new geometry
}


// Return the MBR for a series of geometries
GEOSGeometry *minimumBoundingRectangle(GEOSContextHandle_t ctx,
                                       vector<GEOSGeometry *> geo)
{
    GEOSGeometry *boundingBox = NULL;
    GEOSGeometry *Union, *curMBR;

    for (auto cur = geo.begin(); cur != geo.end(); ++cur)
    {
        curMBR = GEOSEnvelope_r(ctx, *cur);
        if (boundingBox == NULL)
        {
            boundingBox = curMBR;
            continue;
        }

        // Take the union of the current bounding box with the current geometry
        Union = GEOSUnion_r(ctx, boundingBox, curMBR);
        GEOSGeom_destroy_r(ctx,
                           boundingBox); // Have to destroy the old box before
                                         // the new one can be created
        boundingBox = GEOSEnvelope_r(ctx, Union);

        // Free up some memory
        GEOSGeom_destroy_r(ctx, curMBR);
        GEOSGeom_destroy_r(ctx, Union);
    }

    return boundingBox;
}


// Function to calculate the union of local MBRs for a process
GEOSGeometry *calculateAggregateLocalMBR(GEOSContextHandle_t ctx, const vector<GEOSGeometry *> &localMBRs)
{   
    // Combine individual MBRs to get the local MBR for the process
    GEOSGeometry *processLocalMBR = nullptr;

    for (const auto &mbr : localMBRs)
    {       
        if (!processLocalMBR)
        {
            processLocalMBR = GEOSGeom_clone_r(ctx, mbr);
        }  
        else
        { 
            GEOSGeometry *unionResult = GEOSUnion_r(ctx, processLocalMBR, mbr);
            
            GEOSGeom_destroy_r(ctx, processLocalMBR);
            processLocalMBR = unionResult;
        }
         // Free memory for the current MBR
        GEOSGeom_destroy_r(ctx, mbr);
    }

    return processLocalMBR;
}


// Function to calculate the envelope of a geometry
GEOSGeometry *calculateEnvelope(GEOSContextHandle_t ctx, GEOSGeometry *geometry)
{
    return GEOSEnvelope_r(ctx, geometry);
}

// Return the MBR for a single geometry
GEOSGeometry *minimumBoundingRectangle(GEOSContextHandle_t ctx, GEOSGeometry *geometry)
{
    GEOSGeometry *curMBR = GEOSEnvelope_r(ctx, geometry);
    return curMBR;
}

//To create grid over MBR
vector<vector<GEOSGeometry *>> createGrid(GEOSContextHandle_t ctx,
                                          GEOSGeometry *base, int gridSize)
{
    const GEOSCoordSequence *baseCoor =
        GEOSGeom_getCoordSeq_r(ctx, GEOSGetExteriorRing_r(ctx, base));
    unsigned int length;
    GEOSCoordSeq_getSize_r(ctx, baseCoor, &length);

    // Check we got four points
    // Remember that the last one is the same as the first one so that the shape
    // is closed
    if (length != 5)
    {
        aprintf(57 + integerLength(length),
                "Returned Coordinate Sequence has unexpected dimensions! %d\n",
                length);
        return vector<vector<GEOSGeometry *>>();
    }

    double minx, maxx, x, miny, maxy, y;

    // Have to set init values separatly. While -1 would work as
    // a default for the max values, there is not sane default
    // possible for the min values.
    GEOSCoordSeq_getXY_r(ctx, baseCoor, 0, &x, &y);
    minx = x;
    maxx = x;
    miny = y;
    maxy = y;

    for (int i = 1; i < 4; i++)
    {
        // Fetch the values of the current coordinate
        GEOSCoordSeq_getXY_r(ctx, baseCoor, i, &x, &y);

        minx = x < minx ? x : minx;
        maxx = x > maxx ? x : maxx;
        miny = y < miny ? y : miny;
        maxy = y > maxy ? y : maxy;
    }

    double xRange = maxx - minx;
    double yRange = maxy - miny;

    // It feels odd that we have to do this much work to get something as simple
    // as the width and height of a shape.

    double xStep = xRange / gridSize;
    double yStep = yRange / gridSize;

    vector<vector<GEOSGeometry *>> result =
        vector<vector<GEOSGeometry *>>(gridSize);
    for (int r = 0; r < gridSize; r++)
    {
        result[r] = vector<GEOSGeometry *>(gridSize);
        for (int c = 0; c < gridSize; c++)
        {
            result[r][c] = GEOSGeom_createRectangle_r(
                ctx, minx + c * xStep, miny + r * yStep, minx + (c + 1) * xStep,
                miny + (r + 1) * yStep);
        }
    }

    return result;
}

// //to create sketch for polygon within grid
vector<double> sketch(GEOSContextHandle_t ctx,
                      vector<vector<GEOSGeometry *>> *grid, GEOSGeometry *g)

{      

    // The grid should always be square, but just in case it isn't, check the
    // second dimension too
    vector<double> result = vector<double>();

    // cout << "Error 30" << endl;

    int rowNum = grid->size();
    int colNum;
    for (int r = 0; r < rowNum; r++)
    {  //  cout<<"inside rows" <<endl;
        colNum = grid->at(r).size();

        // cout << "Error 30.1" << endl;

        for (int c = 0; c < colNum; c++)
        {
            // cout << "Error 30.1.1" << endl;
            //cout<<"inside columns" <<endl;
            if (grid->at(r).at(c) != NULL)
            {    //cout<<"inside grids" <<endl;
           
                // cout << "Error 30.1.2" << endl;
             double ratio = fillRatio(ctx, grid->at(r).at(c), g);

            //  cout << "Error 30.1.3" << endl;
           //result.push_back(fillRatio(ctx, grid->at(r).at(c), g));
         result.push_back(sqrt(ratio));
           // result.push_back(pow(ratio, 0.5));
                
                
            }
        }

        // cout << "Error 30.2" << endl;

    }

    // cout << "Error 31" << endl;

    result.shrink_to_fit();
    return result;
}





// vector<double> sketch(GEOSContextHandle_t ctx,
//                       vector<vector<GEOSGeometry *>> *grid, GEOSGeometry *g)
// {
//     vector<double> result;
//     int rowNum = grid->size();
//     int colNum;
// 
//     // Initialize min and max values
//     double minVal = std::numeric_limits<double>::max();
//     double maxVal = std::numeric_limits<double>::lowest();
// 
//     // First pass to find min and max
//     for (int r = 0; r < rowNum; r++) {
//         colNum = grid->at(r).size();
//         for (int c = 0; c < colNum; c++) {
//             GEOSGeometry *cell = grid->at(r).at(c);
//             if (cell != NULL) {
//                 double ratio = fillRatio(ctx, cell, g);
//                 if (ratio < minVal) minVal = ratio;
//                 if (ratio > maxVal) maxVal = ratio;
//             }
//         }
//     }
// 
//     // Second pass to calculate normalized values
//     for (int r = 0; r < rowNum; r++) {
//         colNum = grid->at(r).size();
//         for (int c = 0; c < colNum; c++) {
//             GEOSGeometry *cell = grid->at(r).at(c);
//             if (cell != NULL) {
//                 double ratio = fillRatio(ctx, cell, g);
//                 // Apply Min-Max Normalization
//                 double normalizedValue = (ratio - minVal) / (maxVal - minVal);
//                 result.push_back(normalizedValue);
//             }
//         }
//     }
// 
//     result.shrink_to_fit();
//     return result;
// }



// Define the sigmoid scaling function
double sigmoidScale(double value, double k, double x0) {
    // Sigmoid scaling formula
    return 1.0 / (1.0 + exp(-k * (value - x0)));
}


// vector<double> sketch(GEOSContextHandle_t ctx,
//                              vector<vector<GEOSGeometry *>> *grid, GEOSGeometry *g)
// {
//     vector<double> result;
//     int rowNum = grid->size();
//     int colNum;
// 
//     // Parameters for the sigmoid function
//     double k = 10;  // Steepness of the curve
//     double x0 = 0.05;  // Midpoint
// 
//     for (int r = 0; r < rowNum; r++) {
//         colNum = grid->at(r).size();
//         for (int c = 0; c < colNum; c++) {
//             GEOSGeometry *cell = grid->at(r).at(c);
//             if (cell != NULL) {
//                 double ratio = fillRatio(ctx, cell, g);
//                 if (ratio > 0) {
//                     // Apply the sigmoid scaling to the fill ratio only if it's positive
//                     double scaledValue = sigmoidScale(ratio, k, x0);
//                     result.push_back(scaledValue);
//                 } else {
//                     // If ratio is not greater than 0, push back the original ratio
//                     result.push_back(ratio);
//                 }
//             }
//         }
//     }
// 
//     result.shrink_to_fit();
//     return result;
// }


vector<int> LSHHash(vector<double> *sketch, int hashLength,
                    vector<unsigned int> *seeds)

{

// Check for invalid input
    if (!sketch || sketch->empty()) {
        throw std::invalid_argument("LSHHash received an empty sketch vector.");
    }
    if (!seeds || seeds->empty()) {
        throw std::invalid_argument("LSHHash received an empty seeds vector.");
    }
 
    uniform_real_distribution<double> distribution(0.0, static_cast<double>(sketch->size() - 1));
   // uniform_real_distribution<double> distribution(0.0, 0.000001);
    vector<int> hash = vector<int>(hashLength, 0);
   // cout << "LSHHash Debug: sketch size = " << sketch->size() << ", seeds size = " << seeds->size() << endl;


    for (int i = 0; i < hashLength; i++)
    {
        // Seed the random number generator
        minstd_rand generator(seeds->at(i % seeds->size()));
        while (true)
        {
            double r = distribution(generator);
            if (r <= int(r) + (sketch->at(int(r))))
           // cout << "Debug: i = " << i << ", r = " << r << ", r_int = " << r_int << endl;

                break;

            hash[i]++;
        }
    }

    return hash;
}


//generate lsh hash code
vector<int> LSHHashs(vector<double> *sketch, int hashLength, unsigned int seed)
{
    uniform_real_distribution<double> distribution(0.0, static_cast<double>(sketch->size() - 1));
    vector<int> hash = vector<int>(hashLength, 0);

    // Seed the random number generator outside the loop
    minstd_rand generator(seed);

    for (int i = 0; i < hashLength; i++)
    {
        while (true)
        {
            double r = distribution(generator);
            if (r <= int(r) + sketch->at(static_cast<int>(r)))
                break;

            hash[i]++;
        }
       // cout << "Hash[" << i << "] = " << hash[i] << endl;
    }

    return hash;
}



std::vector<std::pair<double, double>> extractCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry)
{
    std::vector<std::pair<double, double>> coordinates;

    int geometryType = GEOSGeomTypeId_r(ctx, geometry);
    if (geometryType != GEOS_POLYGON)
    {
        std::cerr << "Invalid geometry type. Expected Polygon." << std::endl;
        return coordinates;
    }

    // Get the exterior ring of the polygon
    const GEOSGeometry *exteriorRing = GEOSGetExteriorRing_r(ctx, geometry);
    if (exteriorRing == nullptr)
    {
        std::cerr << "Error getting exterior ring." << std::endl;
        return coordinates;
    }

    const GEOSCoordSequence *coordSeq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
    if (coordSeq == nullptr)
    {
        std::cerr << "Error getting coordinate sequence." << std::endl;
        return coordinates;
    }

    unsigned int numPoints;
    GEOSCoordSeq_getSize_r(ctx, coordSeq, &numPoints);

    for (unsigned int i = 0; i < numPoints; ++i)
    {
        double x, y;
        GEOSCoordSeq_getX_r(ctx, coordSeq, i, &x);
        GEOSCoordSeq_getY_r(ctx, coordSeq, i, &y);

        coordinates.push_back(std::make_pair(x, y));
    }

    return coordinates;
    

}


std::pair<std::pair<double, double>, std::pair<double, double>> findBoundingBox(const std::vector<std::pair<double, double>>& points)
{
    if (points.empty()) {
        // Return a default value or handle the case where points are empty
        return {{0.0, 0.0}, {0.0, 0.0}};
    }

    // Initialize the bounding box coordinates
    double xmin = points[0].first;
    double xmax = points[0].first;
    double ymin = points[0].second;
    double ymax = points[0].second;

    // Find the minimum and maximum coordinates
    for (const auto& point : points) {
        xmin = std::min(xmin, point.first);
        xmax = std::max(xmax, point.first);
        ymin = std::min(ymin, point.second);
        ymax = std::max(ymax, point.second);
    }

    // Return the bounding box coordinates
    return {{xmin, ymin}, {xmax, ymax}};
}


GEOSGeometry *createPolygonFromBoundingBox(GEOSContextHandle_t ctx, const std::pair<std::pair<double, double>, std::pair<double, double>> &boundingBox)
{
    // Create a linear ring from the bounding box coordinates
    GEOSCoordSequence *coordSeq = GEOSCoordSeq_create_r(ctx, 5, 2);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 0, boundingBox.first.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 0, boundingBox.first.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 1, boundingBox.second.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 1, boundingBox.first.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 2, boundingBox.second.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 2, boundingBox.second.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 3, boundingBox.first.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 3, boundingBox.second.second);
    GEOSCoordSeq_setX_r(ctx, coordSeq, 4, boundingBox.first.first);
    GEOSCoordSeq_setY_r(ctx, coordSeq, 4, boundingBox.first.second);

    // Create a linear ring and add it to a polygon
    GEOSGeometry *linearRing = GEOSGeom_createLinearRing_r(ctx, coordSeq);
    GEOSGeometry *polygon = GEOSGeom_createPolygon_r(ctx, linearRing, nullptr, 0);

    return polygon;
}

//print coordinates of a geometry
void printCoordinates(GEOSContextHandle_t ctx, const GEOSGeometry *geometry)
{
    int geometryType = GEOSGeomTypeId_r(ctx, geometry);
    if (geometryType != GEOS_POLYGON)
    {
        cerr << "Invalid geometry type. Expected Polygon." << endl;
        return;
    }

    // Get the exterior ring of the polygon
    const GEOSGeometry *exteriorRing = GEOSGetExteriorRing_r(ctx, geometry);
    if (exteriorRing == nullptr)
    {
        cerr << "Error getting exterior ring." << endl;
        return;
    } 
    const GEOSCoordSequence *coordSeq = GEOSGeom_getCoordSeq_r(ctx, exteriorRing);
    if (coordSeq == nullptr)
    {
        cerr << "Error getting coordinate sequence." << endl;
        return;
    }
    unsigned int numPoints;
    GEOSCoordSeq_getSize_r(ctx, coordSeq, &numPoints);

    cout << "Number of Points: " << numPoints << endl;
    for (unsigned int i = 0; i < numPoints; ++i)
    {
        double x, y;
        GEOSCoordSeq_getX_r(ctx, coordSeq, i, &x);
        GEOSCoordSeq_getY_r(ctx, coordSeq, i, &y);

        cout << "Point " << i + 1 << ": (" << x << ", " << y << ")" << endl;
    }
   // cout << "error 7" << endl;

    cout << endl;
}
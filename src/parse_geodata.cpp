// Parse data downloaded from UCR Star

#include "parse_geodata.h"
#include <mpi.h>
#include <unordered_map>


//modified read with wkt for reading both geometry and their unique id
unordered_map<string, GEOSGeometry*> *read_wkt_with_id(const string &filename, GEOSContextHandle_t ctx) {
    auto *geometries = new unordered_map<string, GEOSGeometry*>();
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Unable to open file: " << filename << endl;
        return geometries;
    }

    GEOSWKTReader *reader = GEOSWKTReader_create_r(ctx);

    while (getline(file, line)) {
        size_t firstTabPos = line.find('\t');
        size_t lastTabPos = line.rfind('\t');
        if (firstTabPos == string::npos || firstTabPos == lastTabPos) {
            cerr << "Invalid line format, skipping: " << line << endl;
            continue;
        }

        string id = line.substr(0, firstTabPos);
        string wkt = line.substr(firstTabPos + 1, lastTabPos - firstTabPos - 1);

        GEOSGeometry *geometry = GEOSWKTReader_read_r(ctx, reader, wkt.c_str());
        if (!geometry) {
            cerr << "Invalid WKT for ID " << id << ", skipping" << endl;
            continue;
        }

        // Check if the geometry is a polygon
        if (GEOSGeomTypeId_r(ctx, geometry) != GEOS_POLYGON) {
            cerr << "Geometry for ID " << id << " is not a polygon, skipping" << endl;
            GEOSGeom_destroy_r(ctx, geometry);
            continue;
        }

        // Check the validity of the polygon
        if (GEOSisValid_r(ctx, geometry) != 1) {
            // Attempt to repair the polygon
            GEOSGeometry *validPolygon = GEOSMakeValid_r(ctx, geometry);
            if (validPolygon && GEOSisValid_r(ctx, validPolygon) == 1) {
                GEOSGeom_destroy_r(ctx, geometry);
                geometry = validPolygon;
            } else {
                cerr << "Failed to repair polygon for ID " << id << ", skipping" << endl;
                GEOSGeom_destroy_r(ctx, geometry);
                continue;
            }
        }

        (*geometries)[id] = geometry;
    }

    GEOSWKTReader_destroy_r(ctx, reader);
    file.close();

    return geometries;
}
GEOSGeometry **read_csv(const char *filename, unsigned int *numGeo,
                        GEOSContextHandle_t ctx)
{
    GEOSGeometry **geom = NULL;
    ifstream f(filename);

    if (f.is_open())
    {
        *numGeo = 0;
        unsigned int geoCap = INIT_INPUT_SIZE;
        geom = (GEOSGeometry **)malloc(sizeof(GEOSGeometry *) * geoCap);
        string l;
        GEOSWKTReader *parser = GEOSWKTReader_create_r(ctx);
        // Consumer the header line
        getline(f, l);
        while (getline(f, l))
        {
            // These CSV files are actually tab separated.
            size_t tabIndex = l.find('\t');
            if (tabIndex == string::npos)
            {
                cout << "Cannot find end of first column on " << *numGeo + 1
                     << ", skipping" << endl;
                continue;
            }
            l[tabIndex] = '\0';

            geom[*numGeo] = GEOSWKTReader_read_r(ctx, parser, l.c_str());

            // If the read geometry isn't a POLYGON, drop it.
            if (GEOSGeomTypeId_r(ctx, geom[*numGeo]) != GEOS_POLYGON)
            {
                GEOSGeom_destroy_r(ctx, geom[*numGeo]);
            }
            else
            {
                numGeo++;
            }

            // We may need to double the size of the array
            if (*numGeo >= geoCap)
            {
                geoCap *= 2;
                geom = (GEOSGeometry **)realloc(geom, sizeof(GEOSGeometry *) *
                                                          geoCap);
            }
        }

        GEOSWKTReader_destroy_r(ctx, parser);

        // Trim the length of geom
        geom =
            (GEOSGeometry **)realloc(geom, sizeof(GEOSGeometry *) * (*numGeo));
    }
    else
    {
        cout << "Could not open \"" << filename << "\"" << endl;
    }

    return geom;
}

vector<GEOSGeometry *> *read_wkt(const string filename, GEOSContextHandle_t ctx)
{ 
cout<< filename <<" " << endl;
    vector<GEOSGeometry *> *geom = new vector<GEOSGeometry *>();
    ifstream f(filename);

    if (f.is_open())
    {
        string l;
        GEOSGeometry *possiblePolygon;
        GEOSWKTReader *parser = GEOSWKTReader_create_r(ctx);
        while (getline(f, l))
        {
            // The WKT strings start like '####\tPOLYGON(...)' and we need to
            // start searching on the P
            size_t wktStart = l.find('\t');
            size_t wktEnd = l.rfind('\t');
            if (wktStart == wktEnd)
            {
                aprintf(33 + integerLength(geom->size() + 1),
                        "Cannot find read line %d, skipping\n",
                        geom->size() + 1);
                continue;
            }
            l[wktEnd] = '\0';

            possiblePolygon =
                GEOSWKTReader_read_r(ctx, parser, l.c_str() + wktStart);

            // If the read geometry isn't a POLYGON, drop it.
            if (GEOSGeomTypeId_r(ctx, possiblePolygon) != GEOS_POLYGON)
            {
                GEOSGeom_destroy_r(ctx, possiblePolygon);
            }
            else
            {
                // Check the the polygon is valid. If not, try to repair it and
                // if that fails, destroy it.
                if (GEOSisValid_r(ctx, possiblePolygon) != 1)
                {
                    possiblePolygon = GEOSMakeValid_r(ctx, possiblePolygon);
                    if (GEOSisValid_r(ctx, possiblePolygon) != 1)
                    {
                        aprintf(23, "Failed to fix polygon!\n");
                        GEOSGeom_destroy_r(ctx, possiblePolygon);
                        continue;
                    }
                }
                geom->push_back(possiblePolygon);
            }
        }

        GEOSWKTReader_destroy_r(ctx, parser);

        // Trim the length of geom
        geom->shrink_to_fit();
    }
    else
    {
        aprintf(18 + filename.length(), "Could not open \"%s\"\n",
                filename.c_str());
    }

    return geom;
}

void thread_read_wkt(const string filename, vector<GEOSGeometry *> *geos,
                     mutex *geos_mutex)
{
    // Create a new GEOS context hanlder
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    // Use exisitng read_wkt function
    vector<GEOSGeometry *> *results = read_wkt(filename, ctx);
    aprintf(21 + filename.length() + integerLength(results->size()),
            "Read %d polygons from %s\n", results->size(), filename.c_str());

    // write results
    geos_mutex->lock();

    for (size_t i = 0; i < results->size(); i++)
    {
        geos->push_back(results->at(i));
    }

    geos_mutex->unlock();

    delete results;
    GEOS_finish_r(ctx);
}

vector<GEOSGeometry *> *read_wkt_parallel(const string filename,
                                          unsigned char numThreads,
                                          bool keepTemp)
{
    vector<string> files = splitFile(filename, numThreads);

    vector<thread> *threads = new vector<thread>(numThreads);
    vector<GEOSGeometry *> *geos = new vector<GEOSGeometry *>();
    mutex geos_mutex;
    for (int i = 0; i < numThreads; i++)
    {
        threads->at(i) = thread(thread_read_wkt, files[i], geos, &geos_mutex);
    }

    // Wait for the threads to return and delete them
    for (int i = 0; i < numThreads; i++)
    {
        threads->at(i).join();
        if (!keepTemp)
            unlink(files[i].c_str());
    }
    threads->clear();
    delete threads;
    files.clear();

    return geos;
}

GEOSGeometry **read_geojson(const char *filename, unsigned int *numGeo,
                            GEOSContextHandle_t ctx)
{
    return NULL;
}

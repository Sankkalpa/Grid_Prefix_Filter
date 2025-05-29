#include "thread.h"

void inputThreadFunc(inputIn *in, inputOut *out)
{
    // Create a new GEOS context hanlder
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    out->geos = read_wkt(in->filename, ctx);
    aprintf(21 + in->filename.length() + integerLength(out->geos->size()),
            "Read %d polygons from %s\n", out->geos->size(),
            in->filename.c_str());

    // Continue with data pre-processing, which is centerting all of the
    // geometry and finding the local minimum bounding rectangle.
    for (auto cur = out->geos->begin(); cur != out->geos->end(); ++cur)
        *cur = centerGeometry(ctx, *cur);

    GEOS_finish_r(ctx);
}

// Filter the grid for rarely used cells
gridFilterer gridFilterer::operator()()
{
    int colCount, usage;
    int threadCount = conOut->size();
    int rowCount = grid->size();
    // Iterate over each row
    for (int r = 0; r < rowCount; r++)
    {
        colCount = grid->at(r).size();
        // Iterator over each column
        for (int c = 0; c < colCount; c++)
        {
            // Find total cell usage
            usage = 0;
            for (int i = 0; i < threadCount; i++)
            { 
                usage += conOut->at(i).cellUsage->at(r * colCount + c);
            }

            // Remove the cell
            if (usage < threshold)
            {
                GEOSGeom_destroy_r(ctx, grid->at(r)[c]);
                grid->at(r)[c] = NULL;
                remove->push_back(r * colCount + c);
            }
        }
    }

    int removedCount = remove->size();
    double removedPercent =
        (double)removedCount / (double)(rowCount * colCount);
    aprintf(27 + integerLength(removedCount), "Removed %d (%2.4lf%%) cells!\n",
            removedCount, removedPercent * 100);

    return *this;
}

void constructionThreadFunc(constructionIn *in, constructionOut *out)
{
    // Create a new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    // Generate all of the sketches for the input polygons
    list<vector<double>> sketches = list<vector<double>>();

    for (auto cur = in->geos->begin(); cur != in->geos->end(); cur++)
        sketches.push_front(sketch(ctx, in->grid, *cur));

    // Track cell usage
    out->cellUsage = new vector<int>(sketches.front().size(), 0);

    int idx;
    for (auto sketch = sketches.begin(); sketch != sketches.end(); sketch++)
    {
        idx = 0;
        for (auto element = sketch->begin(); element != sketch->end();
             element++)
        {
            // Count only non-empty /* and non-full */ cells
            if (*element > 0.0 /* && *element < 1.0 */)
                out->cellUsage->at(idx) += 1;
            idx++;
        }
    }

    // Wait for the main thread to process the usage information
    in->bar->arrive_and_wait();

    // Finished wait, revise sketches according to new data
    // While the pointer to removedCells has been here the whole time, the list
    // itself is not populated until after the barrier wait.
    //
    // This also takes advatange that the removedCells list is sorted.
    for (auto sketch = sketches.begin(); sketch != sketches.end(); sketch++)
    {
        int numRemoved = 0;
        int idx = 0;
        auto removed = in->removedCells->begin();
        while (removed != in->removedCells->end())
        {
            if (idx == *removed)
            {
                // "Remove" this element
                numRemoved++;
                removed++;
            }
            else
            {
                // Do not perform the removal, but still copy elements
                sketch->at(idx - numRemoved) = sketch->at(idx);
            }

            idx++;
        }
        sketch->resize(sketch->size() - in->removedCells->size());
    }

    aprintf(28, "Finished trimming sketches!\n");

    // Generate LSH Hashes and insert into the HashMaps
    auto polygon = in->geos->begin();
    for (auto sketch = sketches.begin(); sketch != sketches.end(); sketch++)
    {
        // We need one LSH hash for each table
        for (unsigned int h = 0; h < in->maps->size(); h++)
        {
            vector<int> hash =
                LSHHash(&(*sketch), in->hashLength, &in->seeds->at(h));
            in->maps->at(h)->insert(hash, *polygon);
        }
        polygon++;
    }

    aprintf(29, "Construction Thread Exiting!\n");

    GEOS_finish_r(ctx);
}

void queryThreadFunc(queryIn *in, queryOut *out)
{
    // Create a new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    for (auto cur = in->queryGeos->begin(); cur != in->queryGeos->end(); ++cur)
    {
        // Generate a sketch
        vector<double> s = sketch(ctx, in->grid, *cur);

        // Generate some LSH hashs and query
        vector<pair<double, GEOSGeometry *>> nn;
        for (unsigned int h = 0; h < in->maps->size(); h++)
        {
            vector<int> hash = LSHHash(&s, in->hashLength, &in->seeds->at(h));
            const list<GEOSGeometry *> *bucket = in->maps->at(h)->get(hash);
            if (bucket != NULL)
                for (auto neighbor = bucket->begin(); neighbor != bucket->end();
                     neighbor++)
                {
                    pair<double, GEOSGeometry *> neighborPair = {
                        jaccardDistance(ctx, *cur, *neighbor), *neighbor};
                    // We always have to try to return numNeighbor neighbors
                    if (nn.size() < in->numNeighbor)
                    {
                        nn.push_back(neighborPair);
                        push_heap(nn.begin(), nn.end(), lessthan);
                    }
                    else if (lessthan(neighborPair, nn[0]))
                    {
                        // Found a closer neighbor. One speacial consideration
                        // here is that if the neighbor is closer than the
                        // farthest known nestest neighbor, it is possible that
                        // that neighbor is already in the heap from scanning a
                        // pervious hash table. Remember, close neighbors are
                        // more likely to appear in the same bucket, even across
                        // different hash maps. Since this level of checking is
                        // more expensive, we only do it when their is an actual
                        // chance of adding the shape to the list of nearest
                        // neighbors.
                        //
                        // We can father assume that this is not possible if the
                        // difference from the query object is different. Even
                        // if they are the same, a check of the Jaccard Distance
                        // between the objects is required to check if they are
                        // really the same.

                        bool newNeighbor = true;
                        for (auto n = nn.begin(); n != nn.end(); n++)
                        {
                            if (n->first == neighborPair.first)
                            {
                                if (jaccardDistance(ctx, n->second,
                                                    neighborPair.second) == 0.0)
                                {
                                    newNeighbor = false;
                                }
                            }
                        }

                        if (newNeighbor)
                        {
                            pop_heap(nn.begin(), nn.end(), lessthan);
                            nn[nn.size() - 1] = neighborPair;
                            push_heap(nn.begin(), nn.end(), lessthan);
                        }
                    }
                }
        }
        // Nearest neighbors located, report results
        sort_heap(nn.begin(), nn.end(), lessthan);
        out->results->push_back(
            pair<GEOSGeometry *, vector<pair<double, GEOSGeometry *>>>(*cur,
                                                                       nn));
    }

    GEOS_finish_r(ctx);
}

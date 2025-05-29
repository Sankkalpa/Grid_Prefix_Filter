// utility functions, mostly related to I/O

#ifndef UTILITIES_H
#define UTILITIES_H

#include "geos_c.h"
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <iostream>
#include <list>
#include <mutex>
#include <random>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <unordered_map>


#define OUTPUT_BUFFER_LENGTH 256
#define HASH_VECTOR_MAXIMUM 99

using namespace std;



// ----------------------------------------------------------------------------
// LAMBDA lessthan :
//     Constant anyomous function used to compare elements in the priority
//     queues used in the query threads. Since each element of the queue is a
//     pair, we just need to base the comparison on the first elements of the
//     pairs.
// PARAMETERS :
//     pair<double, GEOSGeometry *> lhs : The left hand side of the "<"
//                                        operation.
//     pair<double, GEOSGeometry *> rhs : The right hand side of the "<"
//                                        operation.
// RETURNS : Boolean, if lhs has a smaller key than rhs.
// ----------------------------------------------------------------------------
const auto lessthan =
    [](pair<double, GEOSGeometry *> lsh, pair<double, GEOSGeometry *> rsh)
{ return lsh.first < rsh.first; };

typedef list<pair<GEOSGeometry *, vector<pair<double, GEOSGeometry *>>>>
    QueryResultSet;

// ----------------------------------------------------------------------------
// FUNCTION aprinf :
//     Write to standard out in an atomic manner. Requires knowing the number
//     of characters to be written since I couldn't find a way to
//     programatically determine the final length of a printf style format
//     string.
// PARAMETERS :
//     int numChar     : The number of characters to print.
//     const char *fmt : printf like format string.
//     ...             : Values for the format specifiers in FMT.
// RETURNS : Nothing
// ----------------------------------------------------------------------------
void aprintf(int numChar, const char *fmt, ...);

// ----------------------------------------------------------------------------
// FUNCTION geosErrorHandler:
//     Write GEOS errors to standard error atomically.
// PARAMETERS :
//     const char *fmt : printf like format string.
//     ...             : Values for the format specifiers in FMT.
// RETURNS : Nothing
// ----------------------------------------------------------------------------
void geosErrorHandler(const char *fmt, ...);

// ----------------------------------------------------------------------------
// FUNCTION :
//     Write GEOS messages to standard out atomically.
// PARAMETERS :
//     const char *fmt : printf like format string.
//     ...             : Values for the format spefiers in FMT.
// RETURNS : Nothing
// ----------------------------------------------------------------------------
void geosMessageHandler(const char *fmt, ...);

// ----------------------------------------------------------------------------
// FUNCTION integerLength :
//     Return the number of characters needed to print an integer.
// PARAMETERS :
//     int i : the integer to be printed.
// RETURNS : The number of characters needed to print an integer.
// ----------------------------------------------------------------------------
int integerLength(int i);

// ----------------------------------------------------------------------------
// FUNCTION splitFile :
//     Use the unix command 'split' to split a file into a set of smaller files
//     in a round-robin style. Used to split large datafiles into smaller files
//     based on where the newlines are.
// PARAMETERS :
//     const string filename  : Relative path to the file to split.
//     unsigned char numSplit : How many smaller files to create.
// RETURNS : List of the names of the smaller files.
// ----------------------------------------------------------------------------
vector<string> splitFile(const string filename, unsigned char numSplit);

// ----------------------------------------------------------------------------
// FUNCTION writerQueryResults :
//     Write the data within RESULTS to files, one file per element in that
//     order. File names will be an integer prepending a type character.
// PARAMETERS :
//     GEOSContextHandle_t ctx : Thread dependent context handle.
//     string filestem         : The location to save the files at.
//     char type               : Letter to append to the integer id in the file
//                               name.
//     QueryResultSet results  : Data to be written.
// ----------------------------------------------------------------------------
void writeQueryResults(GEOSContextHandle_t ctx, string filestem, char type,
                       QueryResultSet results);

// ----------------------------------------------------------------------------
// CLASS HashMap :
//     A one-to-many hash map (meaning one key can have any number of values)
//     which uses a vector of integers as a key. This map has stronger
//     concurrency properties than the standard libraries unordered_map by
//     allowing multiple threads to by able to insert into the map at the same
//     time so long as each is operating on a different bucket.
// MEMBER VARIABLES :
//     void (*elementDeconstructor)
//                           (V value) : Function pointer to a function to
//                                       deconstruct the elements of the map.
//                                       I though I was going to need this,
//                                       but I actually didn't since none of
//                                       the functions in the project operate
//                                       on GEOSGeometry directly, only
//                                       pointers.
//     vector<bucket> table            : Vector used as the basic hash map.
//     vector<int> hashVector          : Randomly generated vector used in the
//                                       hashing process. We take the dot
//                                       product with the key vector then mod
//                                       by the size of the table.
//     int tableSize                   : Length of the TABLE.
//     int currentSize                 : Number of keys in the table.
//     mutex sizeMutex                 : Mutex to protect CURRENTSIZE.
// PRIVATE METHODS :
//     hash         : Generate the hash of a key using a dot product with
//                    HASHVECTOR.
//     increaseSize : Safely increase CURRENTSIZE using SIZEMUTEX.
// METHODS :
//     HashMap()   : Constructor. Takes the following parameters and generates
//                   HASHVECTOR as well as initializing the table.
//                   - int tableSize       : Number of buckets to create.
//                   - int hashSize        : Length of the keys for this map.
//                   - void (*ed)(V value) : element deconstructor
//     ~HashMap () : Deconstructor. Free all of the HashMap memory. Use
//                   ELEMENTDECONSTRUCTOR if given to free the memory of the
//                   elements themselves. However, if the table is storing
//                   pointers this is not needed so long as those pointers
//                   are released elsewhere.
//     size ()     : Return the number of keys in the table.
//     bucketCount : Return the number of buckets in the table.
//     print()     : Print the contents of the table, bucket by bucket, key by
//                   key. The following parameter can modify this behavior.
//                   - bool values : If true print the values, else print only
//                                   the keys in the table.
//     insert()    : Insert a new key-value mapping into the table. If the
//                   required bucket is not avaible, wait for it.
//     tryInsert() : Insert a new key-value mapping into the table. If the
//                   required bucket is not avaible, return false instead.
//     get()       : Return the list of values associated with the input key.
// ----------------------------------------------------------------------------
template <typename V> class HashMap
{
  private:
    typedef struct
    { // Inner level lists (basically vector + list)
        vector<int> *key;
        list<V> *values;
    } keyCollection;
    typedef struct
    { // Top level buckets (basically mutex + list)
        mutex *lock;
        list<keyCollection> *content;
    } bucket;
    void (*elementDeconstructor)(V value);
    vector<bucket> table;
    vector<int> hashVector;
    int tableSize;
    int currentSize;
    mutex sizeMutex;

  //  int hash(vector<int> const key);
  int hash(std::vector<int> const key) const;
    void increaseSize();

  public:
    HashMap(int tableSize, int hashSize, void (*ed)(V value) = nullptr);
    ~HashMap();
    int size();
    int bucketCount();
    void print(bool values);
    void insert(vector<int> const key, V const value);
    bool tryInsert(vector<int> const key, V const value);
    const list<V> *get(vector<int> const key);
    //std::unordered_map<std::vector<int>, int> bucketCountWithKeys(); 
    int getValuesCount(std::vector<int> const &key) const;
    // Method to return the number of buckets in the hash map
    // int size() {
//         std::lock_guard<std::mutex> guard(sizeMutex);
//         return currentSize;
//     }
};

// Possible templates that I will be using. Here so that the compiler creates
// them when compiling this file seperatly from the rest of the project.
template class HashMap<int>;
template class HashMap<GEOSGeometry *>;
template class HashMap<std::string>;
#endif

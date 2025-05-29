#include "util.h"
#include <unordered_map>

// Write to standard error atomically
void geosErrorHandler(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    char buffer[OUTPUT_BUFFER_LENGTH];
    int numChar = vsnprintf(buffer, OUTPUT_BUFFER_LENGTH, fmt, ap);
    buffer[OUTPUT_BUFFER_LENGTH - 2] = '\n';
    write(2, buffer, numChar);
    va_end(ap);
}

// Write to standard out atomically
void geosMessageHandler(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    char buffer[OUTPUT_BUFFER_LENGTH];
    int numChar = vsnprintf(buffer, OUTPUT_BUFFER_LENGTH, fmt, ap);
    buffer[OUTPUT_BUFFER_LENGTH - 2] = '\n';
    write(1, buffer, numChar);
    va_end(ap);
}

// Write to standard out atomically with a known number of characters
void aprintf(int numChar, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    char *buffer = (char *)malloc(sizeof(char) * (numChar + 1));
    vsnprintf(buffer, numChar + 1, fmt, ap);
    write(1, buffer, numChar);
    free(buffer);
    va_end(ap);
}

// Find the number of characters needed to print an integer
int integerLength(int i)
{
    int result;
    if (i <= 0)
    {
        result = 1;
    }
    else
    {
        result = 0;
    }

    i = i < 0 ? -i : i;
    while (i > 0)
    {
        result++;
        i = i / 10;
    }
    return result;
}

vector<string> splitFile(const string filename, unsigned char numSplit)
{
    vector<string> result = vector<string>();
    // Split the file into the correct number of temporary files.
    // Do this using the unix `split` command.
    if (numSplit > 98)
    {
        aprintf(41, "Too many splits requested! numSplit > 98\n");
        return result;
    }

    // If the file is in a different directory, perserve the parent directory
    // structure
    size_t filestem_idx = filename.rfind('/');
    string filestem;
    if (filestem_idx != string::npos)
    {
        // Filestem detected
        filestem = string(filename, 0, filestem_idx + 1);
    }
    else
    {
        // No stem detected, basically the data file is in the current directory
        filestem = string();
    }

    pid_t pid = fork();
    if (pid < 0)
    {
        perror("Unable to fork: ");
        return result;
    }
    else if (pid == 0)
    {
        // Child process
        char chunks[14];
        sprintf(chunks, "--number=r/%u", numSplit);
        execl("/usr/bin/split", "split", chunks, "-d", filename.c_str(),
              filestem.c_str(), (char *)NULL);
    }
    else
    {
        // wait for the split to finish
        wait(NULL);
    }

    // Generate the output of the function
    char fileNumber[3];
    vector<string> files = vector<string>(numSplit);
    for (int i = 0; i < numSplit; i++)
    {
        sprintf(fileNumber, "%02d", i);
        files[i] = string(filestem).append(fileNumber);
    }

    return files;
}



void writeQueryResults(GEOSContextHandle_t ctx, string filestem, char type,
                       QueryResultSet results)
{
    GEOSWKTWriter *writer = GEOSWKTWriter_create_r(ctx);
    char filename[16];
    int idx = 1;
    FILE *f;
    char *polygon;

    for (auto q : results)
    {
        snprintf(filename, 16, "%d%c.txt", idx, type);
        f = fopen(string(filestem).append(filename).c_str(), "w");
        polygon = GEOSWKTWriter_write_r(ctx, writer, q.first);
        fprintf(f, "%s\n", polygon);
        GEOSFree_r(ctx, polygon);

        for (auto n : q.second)
        {
            polygon = GEOSWKTWriter_write_r(ctx, writer, n.second);
            fprintf(f, "%s\n", polygon);
            GEOSFree_r(ctx, polygon);
        }

        fclose(f);
        idx++;
    }

    GEOSWKTWriter_destroy_r(ctx, writer);
}

template <typename V>
HashMap<V>::HashMap(int tableSize, int hashSize, void (*ed)(V value))
    : elementDeconstructor(ed), tableSize(tableSize), currentSize(0)
{
    // Generate values for hashVector
    // Generate values for hashVector
    unsigned int fixedSeed = 123456789; // Choose an appropriate fixed seed value
    minstd_rand generator(fixedSeed); // Initialize minstd_rand with the fixed seed
    uniform_int_distribution<> distrib(0, HASH_VECTOR_MAXIMUM);

    for (int i = 0; i < hashSize; i++)
        hashVector.push_back(distrib(generator));

    // Initiaize table
    table = vector<bucket>();
    for (int i = 0; i < tableSize; i++)
    {
        //		aprintf(21, "Creating new bucket!\n");
        mutex *newMutex = new mutex();
        list<keyCollection> *newList = new list<keyCollection>();
        bucket newBucket = {newMutex, newList};
        table.push_back(newBucket);
    }
}

template <typename V> HashMap<V>::~HashMap()
{
    // Free all allocated memory, moving bucket by bucket
    for (auto b : table)
    {
        b.lock->lock();
        for (auto k : *b.content)
        {
            delete k.key;
            for (auto v : *k.values)
            {
                if (elementDeconstructor != nullptr)
                    elementDeconstructor(v);
            }
            delete k.values;
        }
        delete b.content;
        b.lock->unlock();
        delete b.lock;
    }
}

template <typename V> int HashMap<V>::size()
{
    int result;
    sizeMutex.lock();
    result = currentSize;
    sizeMutex.unlock();
    return result;
}

template <typename V> void HashMap<V>::increaseSize()
{
    sizeMutex.lock();
    currentSize++;
    sizeMutex.unlock();
}

template <typename V> int HashMap<V>::bucketCount() { return tableSize; }


// int HashMap<V>::getValuesCount(vector<int> const &key) const {
//     int idx = hash(key);  // Compute the hash index for the key
//     lock_guard<mutex> guard(*table[idx].lock);  // Lock the bucket
// 
//     for (auto &entry : *table[idx].content) {
//         if (*entry.key == key) {
//             return entry.values->size();  // Return the number of values for the matching key
//         }
//     }
//     return 0;  // Key not found, return 0
// }


template <typename V> void HashMap<V>::print(bool values)
{
    // Print the contens of the entire table
    // Start by locking the whole table so no inserts can be made
    for (auto ptr = table.begin(); ptr != table.end(); ++ptr)
    {
        ptr->lock->lock();
    }

    printf("hashVector = <");
    for (int i : hashVector)
        printf("%2d ", i);
    printf(">\n");

    // Print one bucket at a time
    for (int i = 0; i < tableSize; i++)
    {
        printf("Bucket %03d:\n====================\n", i);
        // Within each bucket their may be multiple keys
        for (auto ptr = table[i].content->begin();
             ptr != table[i].content->end(); ++ptr)
        {
            // Print the key
            // The key itself is a vector
            printf("<");
            for (const int k : *(ptr->key))
                printf("%2d ", k);
            printf(">");

            if (values)
            {
                // Print the values as well, using cout for hopefully automatic
                // type deduction
                printf(": ");
                for (const auto v : *(ptr->values))
                {
                    cout << v << " ";
                }
                printf("\n");
            }
            else
            {
                printf("\n");
            }
        }
        printf("\n");
    }

    // Unlock the table
    for (auto ptr = table.begin(); ptr != table.end(); ++ptr)
    {
        ptr->lock->unlock();
    }
}

template <typename V> int HashMap<V>::hash(vector<int> const key) const
{
    unsigned long long result = 0;
    // char out[256];
    // int len = 11;
    // sprintf(out, "Hashing: ( ");
    for (size_t i = 0; i < hashVector.size(); i++)
    {
        result += key[i] * hashVector[i];
        //	sprintf(out + len, "%d * %d + ", key[i], hashVector[i]);
        //	len += integerLength(key[i]) + integerLength(hashVector[i]) + 6;
    }
    result = result % tableSize;
    // sprintf(out + len - 2, ") %% %d = %d\n", tableSize, result);
    // len += 6 + integerLength(tableSize) + integerLength(result);
    // write(1, out, len);

    return (int)result;
}

template <typename V>
void HashMap<V>::insert(vector<int> const key, V const value)
{
    // Lock the bucket

    int idx = hash(key);
    if (table[idx].lock == NULL)
        aprintf(24 + integerLength(idx), "Illegal mutex value at %d\n", idx);
    table[idx].lock->lock();
    // Scan down the bucket list for matching key
    for (auto ptr = table[idx].content->begin();
         ptr != table[idx].content->end(); ++ptr)
    {
        if (key == *(ptr->key))
        {
            // Insert the new value
            ptr->values->push_front(value);
            table[idx].lock->unlock();
            return;
        }
    }

    // Hash has not been seem before so insert it into the bucket's list
    vector<int> *newKey = new vector(key); // copy constructor, I hope
    list<V> *newValues = new list({value});
    keyCollection newEntry = {newKey, newValues};
    table[idx].content->push_front(newEntry);
    table[idx].lock->unlock();

    increaseSize();
}

template <typename V>
bool HashMap<V>::tryInsert(vector<int> const key, V const value)
{
    // Lock the bucket
    int idx = hash(key);
    if (!table[idx].lock->try_lock())
    {
        return false;
    }

    // Scan down the bucket list for matching key
    for (auto ptr = table[idx].content->begin();
         ptr != table[idx].content->end(); ++ptr)
    {
        if (key == *(ptr->key))
        {
            // Insert the new value
            ptr->values->push_front(value);
            table[idx].lock->unlock();
            return true;
        }
    }

    // Hash has not been seem before so insert it into the bucket's list
    vector<int> *newKey = new vector(key); // copy constructor, I hope
    list<V> *newValues = new list({value});
    keyCollection newEntry = {newKey, newValues};
    table[idx].content->push_front(newEntry);
    table[idx].lock->unlock();

    increaseSize();
    return true;
}

template <typename V> const list<V> *HashMap<V>::get(vector<int> const key)
{
    int idx = hash(key);
    // Scan down the bucket list for matching key
    for (auto ptr = table[idx].content->begin();
         ptr != table[idx].content->end(); ++ptr)
    {
        if (key == *(ptr->key))
        {
            // Return all values associcated with that key
            return ptr->values;
        }
    }

    // key wasn't found
    return nullptr;
}


// // New method to include keys and count the number of polygons associated with each key
// template <typename V>
// std::unordered_map<std::vector<int>, int> HashMap<V>::bucketCountWithKeys() {
//     std::unordered_map<std::vector<int>, int> keyPolygonCount;
// 
//     for (const auto& bucket : table) {
//         std::lock_guard<std::mutex> guard(*bucket.lock);
//         for (const auto& keyColl : *bucket.content) {
//             std::vector<int> key = *keyColl.key;
//             const std::list<V>& polygons = *keyColl.values;
//             keyPolygonCount[key] = polygons.size();
//         }
//     }
// 
//     return keyPolygonCount;
// }

template <typename V>
int HashMap<V>::getValuesCount(std::vector<int> const &key) const {
    int idx = hash(key);  // Compute the hash index for the key
    std::lock_guard<std::mutex> guard(*table[idx].lock);  // Lock the bucket

    for (const auto &entry : *table[idx].content) {
        if (*entry.key == key) {
            return entry.values->size();  // Return the number of values for the matching key
        }
    }
    return 0;  // Key not found, return 0
}

// //to get all the values from bucket
// 
// template <typename V> const list<V> *HashMap<V>::get(vector<int> const key)
//  {
//     int idx = hash(key); // Compute the bucket index based on the hash of the key
//     auto bucket = table[idx].content;
// 
//     if (bucket->empty()) {
//         return nullptr; // If the bucket is empty, return nullptr
//     }
// 
//     // Create a new list to accumulate all values from this bucket
//     list<V>* allValues = new list<V>();
//     for (const auto& entry : *bucket) {
//         allValues->insert(allValues->end(), entry.values->begin(), entry.values->end());
//     }
// 
//     return allValues;
// }
// 
// 


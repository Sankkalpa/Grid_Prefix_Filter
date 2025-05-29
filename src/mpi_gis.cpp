
#include "parse_geodata.h"
#include "query.h"
#include <cmath>
#include <list>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map> 
#include <set>
#include <geos_c.h>
#include "mpi_gis.h"

using namespace std;

vector<string> splitFileMPI(const string filename, MPI_Comm comm) {
    vector<string> result;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // If the file is in a different directory, preserve the parent directory structure
    size_t filestem_idx = filename.rfind('/');
    string filestem;
    if (filestem_idx != string::npos) {
        // Filestem detected
        filestem = string(filename, 0, filestem_idx + 1);
    } else {
        // No stem detected, basically, the data file is in the current directory
        filestem = string();
    }

    // Open the file and get its size
    FILE *file = fopen(filename.c_str(), "rb");
    fseek(file, 0, SEEK_END);
    size_t file_size = ftell(file);
    fclose(file);

    // Calculate the chunk size for each process based on its size
    size_t chunk_size = file_size / size;

    // Calculate the starting position for each process
    size_t start_position = rank * chunk_size;

    // Open the file and read the chunk for the current process
    file = fopen(filename.c_str(), "rb");
    fseek(file, start_position, SEEK_SET);

    // Determine the actual size to read for the last process
    size_t actual_chunk_size = (rank == size - 1) ? (file_size - start_position) : chunk_size;

    // Read the chunk into a buffer
    char *buffer = new char[actual_chunk_size];
    fread(buffer, 1, actual_chunk_size, file);
    fclose(file);

//  // Define the directory where you want to store the temporary files
    string temp_dir = "/home/sp5rt/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/src/data/temp_files";

    // Create a temporary file for each process
    char temp_filename[256];  // Adjust size as needed to accommodate the full path
    snprintf(temp_filename, sizeof(temp_filename), "%stemp_file_%d.txt", temp_dir.c_str(), rank);
    FILE *temp_file = fopen(temp_filename, "wb");
    fwrite(buffer, 1, actual_chunk_size, temp_file);
    fclose(temp_file);

    // 
//     // Create a temporary file for each process
//     char temp_filename[20];
//     snprintf(temp_filename, 20, "temp_file_%d.txt", rank);
//     FILE *temp_file = fopen(temp_filename, "wb");
//     fwrite(buffer, 1, actual_chunk_size, temp_file);
//     fclose(temp_file);

    delete[] buffer;

    // Generate the output of the function
    result.push_back(temp_filename);

    return result;
}
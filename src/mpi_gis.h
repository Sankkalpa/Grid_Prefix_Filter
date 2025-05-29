
#ifndef MPI_GIS_H
#define MPI_GIS_H

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

using namespace std;
vector<string> splitFileMPI(const string filename, MPI_Comm comm);

#endif
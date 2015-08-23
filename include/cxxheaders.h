#ifndef CXXHEADERS_H
#define CXXHEADERS_H

#include "mpi.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using std::endl;
using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::ostringstream;
using std::ios_base;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <algorithm>
using std::max;
using std::min;
using std::swap;
using std::fill_n;

#include <climits>
#include <cfloat>
#include <complex>
using std::complex;
typedef complex<double> dcomplex;

#include <cmath>
using std::abs;
using std::sqrt;
using std::sin;
using std::cos;

#include <cassert>
#include <sys/stat.h>

// Some shortcuts, from Frank Ham
#define FOR_I3 for (int i = 0; i < 3; i++)
#define FOR_J3 for (int j = 0; j < 3; j++)
#define FOR_K3 for (int k = 0; k < 3; k++)
#define FOR_D3 for (int d = 0; d < 3; d++)

#endif

#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <math.h>

#define sgn(x) (double) ((x>0)-(x<0))
#define sqr(x) ((x)*(x))
#define loop(i,a,b) for(int (i)=(a);(i)<(b);(i)++)
#define sz size()
#define ar data()

#define GAMMA_E 0.577215664901532860606512090082l

typedef std::vector<double> vd;

std::ofstream fout;
std::ifstream fin;


#pragma once
using namespace std;

#include <boost/program_options.hpp>

namespace po = boost::program_options;

bool OptionStatus(int argc, char **argv, po::variables_map &VIM);

void ReadVals(po::variables_map &VIM, double &Lx, double &Ly, int &Nx, int &Ny, int &Px, int &Py, double &dt, double &T, double &Re);
#include<cstdlib>
#include<iostream>
using namespace std;

#include <boost/program_options.hpp>

//ProgramOptions.cpp reads input from command line, storing it into variables that are initialised in the functions below. 

namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////
//MEMBER FUNCTIONS
///////////////////////////////////////////////////////////////////

//Checks if the input from command line is in the correct format.
bool OptionStatus(int argc, char **argv, po::variables_map &VIM)
{
   try
   {
    po::options_description desc("Allowed input options");
    // Adding input options with appropriate default values
    desc.add_options()
        ("help", "Produce help message")
        ("Lx", po::value<double>()->default_value(1.0), "Length of the domain in x-direction.")
        ("Ly", po::value<double>()->default_value(1.0), "Length of the domain in y-direction.")
        ("Nx", po::value<int>()   ->default_value(161), "Number of grid points in x-direction.")
        ("Ny", po::value<int>()   ->default_value(161), "Number of grid points in y-direction.")
        ("Px", po::value<int>()   ->default_value(3), "Number of partitions in x-direction (parallel)")
        ("Py", po::value<int>()   ->default_value(2), "Number of partitions in y-direction (parallel)")
        ("dt", po::value<double>()->default_value(0.0005),"Time step size.")
        ("T", po::value<double>() ->default_value(1),"Final time.")
        ("Re", po::value<double>()->default_value(100.0),"Reynolds number.")
    ;
    // Parse command-line arguments and store in buffer VIM
    po::store(po::parse_command_line(argc, argv, desc), VIM);
    po::notify(VIM);
   }
   catch(exception const &e)
   {
    cout << e.what() << endl;
    return false;
   }
   return true;
}

//collects variables from VIM or command line
void ReadVals(po::variables_map &VIM, double &Lx, double &Ly, int &Nx, int &Ny, int &Px, int &Py, double &dt, double &T, double &Re)
{
    Lx = VIM["Lx"].as<double>();
    Ly = VIM["Ly"].as<double>();
    Nx = VIM["Nx"].as<int>();
    Ny = VIM["Ny"].as<int>();
    Px = VIM["Px"].as<int>();
    Py = VIM["Py"].as<int>();
    dt = VIM["dt"].as<double>();
    T =  VIM["T"] .as<double>();
    Re = VIM["Re"].as<double>();
}
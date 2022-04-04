/**
 * Author   :  Priyam Gupta
 * Username :  pg221
 * CID#     :  02110124
 * For      :  AERO70011 High Perfomance Computing Coursework
 *Objective :  Solves Barkley Model using central finite differences for spatial discritization
               and forward Euler for Time Intergration.
 **/


#include <iostream>                         // For Standard input/output
using namespace std;
#include <iomanip>                          // Matrix Output Spacing
#include <fstream>                          // File Handling
#include <omp.h>                            // For parallel run
#include "ReactionDiffusion.h"              // ReactionDiffusion Class
#include <string>                           // String input for File Name
#include <boost/program_options.hpp>        // Boost Library for terminal input
namespace po = boost::program_options;





int main(int argc, char* argv[]){

    //Boost program_options
    po::options_description opts("Available options.");
    opts.add_options()
    ("dt"       ,po::value<double>()->default_value(0.001)    ,"Time-step to use."         )
    ("T"        ,po::value<double>()->default_value(100)      ,"Total Integration Time."   )
    ("Nx"       ,po::value<int>()   ->default_value(101)      ,"Number of grid points in x")
    ("Ny"       ,po::value<int>()   ->default_value(101)      ,"Number of grid points in y")
    ("a"        ,po::value<double>()->default_value(0.75)     ,"Value of parameter a"      )
    ("b"        ,po::value<double>()->default_value(0.06)     ,"Value of parameter b"      )
    ("mu1"      ,po::value<double>()->default_value(5.0)      ,"Value of parameter mu1"    )
    ("mu2"      ,po::value<double>()->default_value(0.0)      ,"Value of parameter mu2"    )
    ("eps"      ,po::value<double>()->default_value(50.0)     ,"Value of parameter eps"    )
    ("np"       ,po::value<double>()->default_value(1)      ,"Number of threads to be spawned"    )
    ("filename" ,po::value<string>()->default_value("uv.txt") ,"Name of Written uv File."  )
    ("help", "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);
    if (vm.count("help")) {
      std::cout << opts << std::endl;
      return 0;
    }

    const double dt  = vm["dt"].as<double>();
    const double T   = vm["T"].as<double>();
    const int    Nx  = vm["Nx"].as<int>();
    const int    Ny  = vm["Ny"].as<int>();
    const double a   = vm["a"].as<double>();
    const double b   = vm["b"].as<double>();
    const double mu1 = vm["mu1"].as<double>();
    const double mu2 = vm["mu2"].as<double>();
    const double eps = vm["eps"].as<double>();
    const int    nthreads = vm["np"].as<double>();
    const string filename = vm["filename"].as<string>();
    const double dx  = 1;
    const double dy  = 1;


    //Instantiating an Object of the ReactionDiffusion Class
    ReactionDiffusion rd;

    //Setting parameters
    rd.SetParameters(Nx, Ny, a, b, eps, mu1, mu2, T, dt, dx, dy, nthreads);

    //Setting Initial Conditions in U and V
    rd.SetInitialConditions();

    //Integrating the PDE over Time
    rd.TimeIntegrate();

    //Writing U and V values into a file
    rd.write(filename);

}

/*
 * This file is part of the Quasi-Static Time-Dynamic Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <string>
#include "../src/headers/math33.H"
#include "../src/headers/single.H"
#include "../src/headers/cup.H"

/*
g++ -Wall -I/usr/include/ -L/usr/local/lib time_behavior.cxx -o tim -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main(int argc, char** argv) {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
    double   omeeV, omemi, omema, E0, *fro, tpump, T, eps_b;
    char mtl[16], mdl[16], sol[16], active[16];
    if (argv[1]==0){
        cout<<endl<<"  Usage: "<<argv[0]<<" <omega in eV>"<<endl<<endl;
        exit(0);
        }
    omeeV=atof(argv[1]);
    nanosphere  simulation;
    simulation.init();
    
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        std::cerr << "Error: Cannot open input file" << std::endl;
        return 1;
    }
    ifstream time("../data/input/time.dat");
    if (!time) {
        std::cerr << "Error: Cannot open input file" << std::endl;
        return 1;
    }
    
    nano>>simulation.r1>>simulation.Dome>>simulation.ome_0>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0;
    time>>T>>tpump;  
        
    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps_b=simulation.set_host(sol);
    
    // Perform the time_behavior calculation

 
    fro=simulation.frohlich_optimal(omemi, omema, eps_b);
    
    // Inform the user about the test
    cout << "Calculating the time_behavior up to "<<T<<" ps\n";
    cout << "                switching the pump on at "<<tpump<<"ps ...\n\n";
    cout << "* Beware if  "<<simulation.G<<" is greater than "<<fro[1]<<"\n";
    cout << "* you can get exponential behaviors in the\n";
    cout << "* analytical results for some frequencies\n\n";

    cout << "Parameters:\n";
    cout << "  Metal model: " << mdl << "\n";
    cout << "  Metal type: " << mtl << "\n";
    cout << "  Solvent: " << sol << "\n\n";

    cout << "Running analytical calculation...\n";
    simulation.analytical(mdl, mtl, sol, E0, omeeV, T, tpump);
    cout << "Running numerical calculation...\n";
    simulation.numerical(mdl, mtl, sol, E0, omeeV, T, tpump);

    // Output the results

    // Save the results to a file
    ofstream output("results/time_behavior.log");
    if (output.is_open()) {
        output << "Steady-state polarizability calculation results:\n";
	output << "  Frequency: " << omeeV << " eV\n";
	output << "  Metal model: " << mdl << "\n";
        output << "  Metal type: " << mtl << "\n";
        output << "  Time range: [0, " << T << "] ps\n";
        output << "  The pump is switched on at "<<tpump<<" ps\n";
        output << "  Solvent: " << sol << "\n";
        output << "\nCheck simulation output for detailed results.\n";
        output.close();
        cout.precision(10);
        cout.setf(ios::fixed);
        cout <<"> Output saved in the following files:"<<endl;
	    cout<<">   ../data/output/anltime.dat"<<endl;
	    cout<<">   ../data/output/anlfunc.dat"<<endl;
	    cout<<">   ../data/output/numtime.dat"<<endl;
	    cout<<">   ../data/output/numfunc.dat"<<endl;
        cout<<"> log saved in results/time_behavior.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}

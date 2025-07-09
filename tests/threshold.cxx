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
g++ -Wall -I/usr/include/ -L/usr/local/lib threshold.cxx -o trs -lgsl -lgslcblas -lm -larmadillo
*/

using namespace std;

int main() {
    // Initialize variables for input parameters
    // Create an instance of the nanosphere class
    double   omemi, omema, E0, *fro, eps_b;
    char mtl[16], mdl[16], sol[16], active[16];
    
    nanosphere  simulation;
    simulation.init();
    
    ifstream nano("../data/input/nanosphere_eV.dat");
    if (!nano) {
        std::cerr << "Error: Cannot open input file" << std::endl;
        return 1;
    }

    nano>>simulation.r1>>simulation.Dome>>simulation.ome_0>>simulation.G>>omemi>>omema>>mtl>>mdl>>active>>sol>>E0;
    
    // Inform the user about the test
    cout << "Calculating the threshold gain *G_th* and the\n";
    cout << "                threshold frequency *ome_th*...\n";
    cout << "Parameters:\n";
    cout << "  Metal model: " << mdl << "\n";
    cout << "  Metal type: " << mtl << "\n";
    cout << "  Spectral range: [" << omemi << ", " << omema << "] eV\n";
    cout << "  Solvent: " << sol << "\n\n";

    simulation.set_metal(mtl,mdl,1);
    simulation.set_active(active);
    eps_b=simulation.set_host(sol);
    
    // Perform the threshold calculation
    cout << "Calling frolich subroutine...\n";
 
    fro=simulation.frohlich_optimal(omemi, omema, eps_b);

    simulation.steady_state(mdl, mtl, sol, omemi, omema, 10000);

    // Output the results

    // Save the results to a file
    ofstream output("results/threshold.log");
    if (output.is_open()) {
        output << "Steady-state polarizability calculation results:\n";
        output << "  Metal model: " << mdl << "\n";
        output << "  Metal type: " << mtl << "\n";
        output << "  Spectral range: [" << omemi << ", " << omema << "] eV\n";
        output << "  Solvent: " << sol << "\n";
        output << "\nCheck simulation output for detailed results.\n";
        output.close();
        cout.precision(10);
        cout.setf(ios::fixed);
        cout <<"> G_th   : "<<fro[1]<<endl;
	    cout <<"> ome_th : "<<fro[0]<<" eV"<<endl;
        cout<<"> log saved in results/threshold.log"<<endl;
    } else {
        cerr << "Error: Could not open file for writing results.\n";
    }

    return 0;
}

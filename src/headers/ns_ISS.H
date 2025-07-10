/*
 * This file is part of the Nano-Shell Simulation Project.
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

std::complex<double> img = std::complex<double> (0., 1.);
    
void complete(std::vector<std::pair<double, double>>& vec, double omemi, double omema, double dome, double value=0.) {
    vec.insert(vec.begin(), std::make_pair(vec[0].first - dome, value));
    vec.insert(vec.begin(), std::make_pair(omemi, value));

    vec.push_back(std::make_pair(vec[vec.size() - 1].first + dome, value));

    vec.push_back(std::make_pair(omema, value));
    }
        
std::vector<std::pair<double,std::complex<double>>> gimme_emi_kap(nanosphere ns, char* mdl, char* mtl, char* hst, double omemi, double omema, int omeN=1000, char* sol=0, double rho=0){
        std::vector<std::pair<double,std::complex<double>>> vkamp;
        int omi;
        std::complex<double> GG, OmeH, GamP, OmeP, **A, *p0, *p1, *p2, gam1, gam2, xi, *kap, *odv;
        
        double *nv, omeeV, dome  = (omema-omemi)/omeN, ome;
        double  eps_b, eps_s, tau2, gamd, ome0;

        if ((sol != NULL)) {
            eps_s=ns.set_host(sol);
            } else eps_s=0;
        eps_b=ns.set_host(hst);

        A = new std::complex<double>*[3];
        for(int i = 0; i <= 2; i++)
            A[i] = new std::complex<double>[3];
        
        kap  = new std::complex<double>[3];
        nv   = new double[4];
        odv  = new std::complex<double> [3];

        // normalized variables
        nv = ns.normalized_variables ();
        tau2 = nv[1];
        gamd = nv[2];
        ome0 = nv[3];
        GG   = ns.set_GG(ns.G, tau2);

        p0 = pcfc0(ns.ceps_inf, eps_b, eps_s, rho);
        p1 = pcfc1(ns.ceps_inf, eps_b, eps_s, rho);
        p2 = pcfc2(ns.ceps_inf, eps_b, eps_s, rho);
        
        for (omi=0; omi<=omeN; omi++){
            omeeV = omemi + omi*dome;
        
            ome  = omeeV/ns.Ome_p;

            odv  = ns.set_ome_dep_vrbls(ome, ome0, tau2, gamd);
            OmeH = odv[0];
            OmeP = odv[1];
            GamP = odv[2];
    
            A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);
            kap  = eigenvalues(A);
            
            double maxVal = kap[0].real();
            int    imax = 0;

            for (int i = 1; i < 3; i++) {
                if (kap[i].real() > maxVal) {
                    maxVal = kap[i].real();
                    imax = i;
                    }
                }
            vkamp.push_back(std::make_pair(omeeV, kap[imax]));
            }
        return vkamp;
        }



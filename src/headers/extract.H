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

std::complex<double> img = std::complex<double> (0., 1.);
    
void complete(std::vector<std::pair<double, double>>& vec, double omemi, double omema, double dome, double value=0.) {
    vec.insert(vec.begin(), std::make_pair(vec[0].first - dome, value));
    vec.insert(vec.begin(), std::make_pair(omemi, value));

    vec.push_back(std::make_pair(vec[vec.size() - 1].first + dome, value));

    vec.push_back(std::make_pair(omema, value));
    }
        

double interpolate(double x1, double x2, double f1, double f2, double x) {
    // Perform linear interpolation
    return f1 + (f2 - f1) * (x - x1) / (x2 - x1);
}

std::pair<double, double> find_zeros(const std::vector<double>& x, const std::vector<double>& fx, double epsilon = 1e-6, int max_iterations = 100) {
    double zero_1 = 0.0;
    double zero_2 = 0.0;
    double c, fc;
    
    for (int i = 0; i < int(x.size()); i++) {
        double x1 = x[i];
        double x2 = x[i + 1];
        double f1 = fx[i];
        double f2 = fx[i + 1];
        
        
        if (f1 * f2 < 0) {
            double a = x1;
            double b = x2;
            int iterations = 0;
            while (std::abs(b - a) > epsilon && iterations < max_iterations) {
                iterations++;
                c = (a + b) / 2.0;
                fc = interpolate(x1, x2, f1, f2, c);

                if (f1 * fc < 0) {
                    b = c;
                }
                else {
                    a = c;
                }
            }

            if (zero_1 == 0.0) {
                zero_1 = c;
            }
            else {
                zero_2 = c;
                break;
            }
        }
        if (i==int(x.size()) - 1 && zero_2 == 0.0) {
            zero_2 = x[i];
            break;
            }
        
    }
    return std::make_pair(zero_1, zero_2);
}

std::vector<double> extract_ome(std::vector<std::pair<double,std::complex<double>>> vectorOfPairs){
    std::vector<double> extractFirst;
        int i =0;
        for (const auto& pair : vectorOfPairs) {
            i++;
            extractFirst.push_back(pair.first);
            if (i==int(vectorOfPairs.size())-1) break;
        }
        return extractFirst;
        }
    
std::vector<double> extract_ralph(std::vector<std::pair<double,std::complex<double>>> vectorOfPairs){
    std::vector<double> extractSecnd;
    int i =0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractSecnd.push_back(real(pair.second));
        if (i==int(vectorOfPairs.size())-1) break;
        }
    return extractSecnd;
    }

std::pair<double, double> fnd_extrm(std::vector<std::pair<double,std::complex<double>>> vkape, double Ome_p){
    std::vector<double> tmp;
    double extr_1, extr_2;
    int never=1;
    for (int ii=0; ii<int(vkape.size()); ii++){
        if (real(vkape[ii].second)>=0) {
            tmp.push_back(vkape[ii].first-imag(vkape[ii].second)*Ome_p);
            never=0;
            }
        }
    sort(tmp.begin(), tmp.end());
    if (never ==0 ){
        extr_1=tmp[0];
        extr_2=tmp[tmp.size()-1];
        } else{
        extr_1=666;
        extr_2=777;
        } 
    return std::make_pair(extr_1, extr_2);
    }
    
std::vector<double> extract_ialph(std::vector<std::pair<double,std::complex<double>>> vectorOfPairs){
    std::vector<double> extractSecnd;
    int i =0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractSecnd.push_back(imag(pair.second));
        if (i==int(vectorOfPairs.size())-1) break;
        }
    return extractSecnd;
    }

bool compareSecond(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return fabs(a.second) > fabs(b.second);
}

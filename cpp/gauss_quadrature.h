//
//  gauss_quadrature.h
//  
//
//  Created by Jacob Koester on 1/5/15.
//
//

#ifndef ____gauss_quadrature__
#define ____gauss_quadrature__

#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> gauss_pts(const int order);
std::vector<double> gauss_wts(const int order);
std::tuple<std::vector<double>, std::vector<double>> pts_and_wts(const int order, double *interval_ends, bool debug = false);
std::tuple<std::vector<double>, std::vector<double>> pts_and_wts(const int order, bool debug = false);
double integrate(const std::vector<double> & values, const std::vector<double> & weights);
void test();

#endif /* defined(____gauss_quadrature__) */

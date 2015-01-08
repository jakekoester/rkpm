//
// quadrature_driver.cpp
//
// Created by Jacob Koester on 1/5/15
//
// Simple program to test function in qauss_quadrature.cpp
//

#include "quadrature_driver.h"

int main(int argc, char* argv[])
{
    // Test some of the functions in gauss_quadrature.cpp
    
    std::tuple<std::vector<double>, std::vector<double>> p_w;
    
    double interval_pts[6] = {-2., 2.};
    p_w = pts_and_wts(3,interval_pts,true);
    p_w = pts_and_wts(3,true);
    
    test();
}
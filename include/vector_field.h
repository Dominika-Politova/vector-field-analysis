#pragma once
#include <vector>
#include <iostream>
#include <cmath>

struct vector_3d {
    double x, y, z;

    vector_3d operator + (const vector_3d& v);
    vector_3d operator - (const vector_3d& v);
    
    friend std::ostream& operator << (std::ostream& os, const vector_3d& v);
};

class grid {
private:
    double dx, dy, dz; // Шаг сетки
    int nx, ny, nz; // Количество элементов на сетке
    std::vector<vector_3d> data;

    static void thomas_algorithm(const std::vector<double>& a,
                                 const std::vector<double>& b,
                                 const std::vector<double>& c,
                                 const std::vector<double>& d,
                                 std::vector<double>& x);

    static std::vector<double> poisson_equation_integrate(const std::vector<double>& div_slice, double h);

public:
    grid(int nx = 5, int ny = 5, int nz = 5, 
         double dx = 0.1, double dy = 0.1, double dz = 0.1);
    
    void get_coords(int i, int j, int k, double& x, double& y, double& z);
    
    vector_3d& operator()(int i, int j, int k);
    
    double div(int i, int j, int k);
    
    vector_3d rot(int i, int j, int k);
    
    void analyze_field();
};
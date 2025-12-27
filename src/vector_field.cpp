#include "../include/vector_field.h"
#include <iostream>
#include <cmath>

using namespace std;

// Реализация методов vector_3d
vector_3d vector_3d::operator + (const vector_3d& v) {
    return {x + v.x, y + v.y, z + v.z};
}

vector_3d vector_3d::operator - (const vector_3d& v) {
    return {x - v.x, y - v.y, z - v.z};
}

ostream& operator << (ostream& os, const vector_3d& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

// Реализация методов grid

// Интегрирование уравнения Пуассона (для одномерного случая) для нахождения потенциала
vector<double> grid::poisson_equation_integrate(const vector<double>& div_slice, double h) {
    int n = div_slice.size(); // div_slice - срез дивергенции вдоль оси x
    vector<double> phi(n, 0.0);
    vector<double> F_pot_temp(n, 0.0);

    for (int i = 1; i < n; i++) {
        double avg_div = 0.5 * (div_slice[i - 1] + div_slice[i]); // Метод трапеций для интегрирования
        F_pot_temp[i] = F_pot_temp[i - 1] + avg_div * h;
    }

    for (int i = 1; i < n; i++) {
        double avg_F = 0.5 * (F_pot_temp[i - 1] + F_pot_temp[i]);
        phi[i] = phi[i - 1] - avg_F * h;
    }

    return phi;
}

grid::grid(int nx, int ny, int nz, double dx, double dy, double dz) {
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    data.resize(nx * ny * nz);
}

void grid::get_coords(int i, int j, int k, double& x, double& y, double& z) {
    x = i * dx;
    y = j * dy;
    z = k * dz;
}

vector_3d& grid::operator()(int i, int j, int k) {
    return data[i * ny * nz + j * nz + k]; // Используется линейная индексация
}

// вычисление дивергенции и ротора векторного поля
// для вычисления производных используются центральные разности второго порядка
// на границах возвращается 0
double grid::div(int i, int j, int k) {
    if (i <= 0 || i >= nx - 1 || j <= 0 || j >= ny - 1 || k <= 0 || k >= nz - 1) {
        return 0.0;
    }
    double dP_x = (operator()(i + 1, j, k).x - operator()(i - 1, j, k).x) / (2 * dx);
    double dQ_y = (operator()(i, j + 1, k).y - operator()(i, j - 1, k).y) / (2 * dy);
    double dR_z = (operator()(i, j, k + 1).z - operator()(i, j, k - 1).z) / (2 * dz);
    return (dP_x + dQ_y + dR_z);
}

vector_3d grid::rot(int i, int j, int k) {
    vector_3d rotor = {0.0, 0.0, 0.0};
    if (i <= 0 || i >= nx - 1 || j <= 0 || j >= ny - 1 || k <= 0 || k >= nz - 1) {
        return rotor;
    }
    double dP_y = (operator()(i, j + 1, k).x - operator()(i, j - 1, k).x) / (2 * dy);
    double dP_z = (operator()(i, j, k + 1).x - operator()(i, j, k - 1).x) / (2 * dz);
    double dQ_x = (operator()(i + 1, j, k).y - operator()(i - 1, j, k).y) / (2 * dx);
    double dQ_z = (operator()(i, j, k + 1).y - operator()(i, j, k - 1).y) / (2 * dz);
    double dR_x = (operator()(i + 1, j, k).z - operator()(i - 1, j, k).z) / (2 * dx);
    double dR_y = (operator()(i, j + 1, k).z - operator()(i, j - 1, k).z) / (2 * dy);

    rotor.x = (dR_y - dQ_z);
    rotor.y = (dP_z - dR_x);
    rotor.z = (dQ_x - dP_y);

    return rotor;
}

// Анализ векторного поля с разложением Гельмгольца: F = F_potential + F_solenoidal
// Метод анализирует поле вдоль центральной линии по оси x (j = ny/2, k = nz/2)
void grid::analyze_field() {
    int center_j = ny / 2;
    int center_k = nz / 2;

    vector<double> div_slice(nx, 0.0);
    for (int i = 1; i < nx - 1; i++) {
        div_slice[i] = (operator()(i + 1, center_j, center_k).x -
                       operator()(i - 1, center_j, center_k).x) / (2 * dx);
    }

    // на границах используются односторонние разности
    div_slice[0] = (operator()(1, center_j, center_k).x - operator()(0, center_j, center_k).x) / dx;
    div_slice[nx - 1] = (operator()(nx - 1, center_j, center_k).x - operator()(nx - 2, center_j, center_k).x) / dx;

    vector<double> phi = poisson_equation_integrate(div_slice, dx);
    vector<vector_3d> F_potential(nx);

    for (int i = 1; i < nx - 1; i++) {
        F_potential[i].x = -(phi[i + 1] - phi[i - 1]) / (2 * dx); //F_potential = -grad ф
        F_potential[i].y = 0.0;
        F_potential[i].z = 0.0;
    }

    F_potential[0].x = -(phi[1] - phi[0]) / dx;
    F_potential[0].y = 0.0;
    F_potential[0].z = 0.0;

    F_potential[nx - 1].x = -(phi[nx - 1] - phi[nx - 2]) / dx;
    F_potential[nx - 1].y = 0.0;
    F_potential[nx - 1].z = 0.0;


    // F_solenoidal = F - F_potential
    vector<vector_3d> F_solenoidal(nx);
    for (int i = 1; i < nx - 1; i++) {
        F_solenoidal[i] = operator()(i, center_j, center_k) - F_potential[i];
    }

     cout << " (Only x-coordinate!) " << endl;
    cout << "\nchange of the rotor at several points" << endl;
    for (int i = 4; i < 7; i++) {
        vector_3d rot_potential = rot(i, center_j, center_k);
        cout << "  i=" << i << ": rot = " << rot_potential << endl;
    }

    double div_sol = (F_solenoidal[nx / 2 + 1].x -
                     F_solenoidal[nx / 2 - 1].x) / (2 * dx);
    cout << "div(F_sol) in center:     " << div_sol << " (must be 0)" << endl;

    cout << "\nresults in center (i=" << nx / 2 << ", x=" << (nx / 2) * dx << "):" << endl;
    cout << "  full field: " << operator()(nx / 2, center_j, center_k) << endl;
    cout << "  potential part: " << F_potential[nx / 2] << endl;
    cout << "  solenoidal part: " << F_solenoidal[nx / 2] << endl;
    cout << "  sum of the solenoidal and potential parts: " << F_potential[nx / 2] + F_solenoidal[nx / 2] << endl;
    cout << "  error: "
         << operator()(nx / 2, center_j, center_k) - (F_potential[nx / 2] + F_solenoidal[nx / 2]) << endl;
    cout << "   div(F) = " << div_slice[nx / 2] << endl;
    cout << "   rot(F) = " << rot(nx / 2, center_j, center_k) << endl;
}

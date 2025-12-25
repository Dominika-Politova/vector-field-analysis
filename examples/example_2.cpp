#include "../include/vector_field.h"
#include <iostream>
#include <cmath>

using namespace std;

int main() {
    int n = 100;
    double h = 0.01;
    grid Grid(n, n, n, h, h, h);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                double x, y, z;
                Grid.get_coords(i, j, k, x, y, z);
                Grid(i, j, k) = {x * x, 0, 0}; // 100% potential field
            }
        }
    }

    cout << "Divergence at (10,10,10): " << Grid.div(10, 10, 10) << endl;

    vector_3d rotation = Grid.rot(10, 10, 10);
    cout << "Rotation at (10,10,10): (" 
         << rotation.x << ", " 
         << rotation.y << ", " 
         << rotation.z << ")" << endl;

    cout << "\n    Field Analysis    " << endl;
    Grid.analyze_field();

    return 0;
}
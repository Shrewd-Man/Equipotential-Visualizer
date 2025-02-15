//
//  main.c
//  Equipotential-Visualizer
//
//  Created by Kenneth Anderson on 2/15/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// define the spacial boundries of the simulation
#define Xmin -5.0
#define Xmax 5.0
#define Ymin -5.0
#define Ymax 5.0

// define the resolution (points) of the simulation within the boundries
#define Nx 100
#define Ny 100

// define the placement and magnitude of the point charge
#define chargeX 0
#define chargeY 0
#define chargeMag (1.062e-19)

// define coulombs constant
#define coul (8.99e9)

double V[Nx][Ny];

/**
 *  @brief Function calculates the potential difference of any provided point
 *
 *  The function uses the universal formula and currently computes that based on the points distance to (0,0), assuming theres a charge with e at (0,0)
 *
 *  @param Px is the given point x value of the place to compute the difference
 *  @param Py is the given point y value of the place to compute the difference
 *
 *  @return Func rerurns the value of the potential difference at the given point
 *
 *  @note The function currently assumes there is only one charge present at (0,0) with magnitude e
 */
double calcVolt(int Px, int Py) {
    
    double dx = (Xmax - Xmin) / (Nx - 1);
    double dy = (Ymax - Ymin) / (Ny - 1);
    
    double xGridPoint = Xmin + Px * dx;
    double yGridPoint = Ymin + Py * dy;
    
    double distanceToCharge = sqrt(pow(fabs(xGridPoint - chargeX), 2) + pow(fabs(yGridPoint - chargeY), 2));
    
    return (coul * chargeMag)/distanceToCharge;
}

void setupGrid(void) {
    double dx = (Xmax - Xmin) / (Nx - 1);
    double dy = (Ymax - Ymin) / (Ny - 1);

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            double x = Xmin + i * dx;
            double y = Ymin + j * dy;
            V[i][j] = calcVolt(i, j);
            printf("Grid Point (%d, %d) -> (%.2f, %.2f) with ∆V of %.10e\n", i, j, x, y, V[i][j]);
        }
    }
}

int main(int argc, const char * argv[]) {
    setupGrid();
    
    return 0;
}

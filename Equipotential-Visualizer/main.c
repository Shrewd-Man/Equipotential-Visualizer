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
#define Nx 1000
#define Ny 1000

// define the placement and magnitude of the point charge
#define chargeX 0
#define chargeY 0
#define chargeMag (1.062e-19)

// define coulombs constant
#define coul (8.99e9)

// Initialize 2 dimensional array and the "dx & dy" values.
double V[Nx][Ny];
double dx = (Xmax - Xmin) / (Nx - 1);
double dy = (Ymax - Ymin) / (Ny - 1);

/**
 *  @brief Function calculates the potential difference of any provided point
 *
 *  The function uses the universal formula and currently computes that based on the points distance to (0,0), assuming theres a charge with e at (0,0)
 *
 *  @param Px is the given point x value of the place to compute the difference
 *  @param Py is the given point y value of the place to compute the difference
 *
 *  @return Func returns the value of the potential difference at the given point
 *
 *  @note The function currently assumes there is only one charge present at (0,0) with magnitude e
 */
double calcVolt(int Px, int Py) {
    
    double xGridPoint = Xmin + Px * dx;
    double yGridPoint = Ymin + Py * dy;
    
    double distanceToCharge = sqrt(pow(fabs(xGridPoint - chargeX), 2) + pow(fabs(yGridPoint - chargeY), 2));
    
    //printf("Grid Point (%d, %d) -> (%.2f, %.2f)\n", Px, Py, xGridPoint, yGridPoint); <-- Save for printing in case of error
    
    return (coul * chargeMag)/distanceToCharge;
}

// The setupGrid function loops through each element in the 2d array and calculates the potential difference of the point
void setupGrid(void) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            V[i][j] = calcVolt(i, j); // Use to calculate the difference
        }
    }
}

int findFurthestWall(int Cx, int Cy, int indexReq) { // find the wall furthest from the charge
    int distances[4];
    int toRightWall = (Nx - Cx);
    int toLeftWall = Cx;
    int toTopWall = Cy;
    int toBottomWall = (Ny - Cy);
    
    distances[0] = toTopWall;
    distances[1] = toBottomWall;
    distances[2] = toRightWall;
    distances[3] = toLeftWall;
    
    int maxIndex = 1;
    for (int i = 2; i <= 4; i++) {
        if (distances[i] > distances[maxIndex]) {
            maxIndex = i;
        }
    }
    
    // return value based on indexReq
    if(indexReq) {
        return maxIndex;
    } else {
        return distances[maxIndex];
    }
}

/**
 * @note The plan is to use an array thats sized based on the distance from the charge to the furthest wall.
 * The function will then look through each point on the line from the charge to the furthest wall and take the value of every 25 points, storing it in the array.
 * Using each value, the program will then render the points where the values fall within a specified "grace" range of one of the values in the array as a gray color,
 * and the array index the value falls under will decide the shade of gray the values should be.
 *
 * Later, use a for loop to iterate this for each charge using the same array, that way we dont need to allocate a ton of memory, and we can just use the one array to calculate lines for every charge
 */

// render
void renderEqLines(int Cx, int Cy) {
    printf("Pass\n");
    
    // Use 0 for up, 1 for down, 2 for right, and 3 for left
    int directionToFurthestWall = findFurthestWall(Cx, Cy, 1); // gather the direction of the distance
    int distanceToWall = findFurthestWall(Cx, Cy, 0); // gather the VALUE of the distance
    
    int n = distanceToWall / 25;
    double *lineData = (double *)malloc(n * sizeof(double));
    
    free(lineData);
}

int main(int argc, const char * argv[]) {
    int Cx, Cy;
    
    Cx = 50;
    Cy = 50;
    
    setupGrid();
    
    renderEqLines(Cx, Cy);
    
    return 0;
}

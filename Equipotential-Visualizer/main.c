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

// define the resolution (points) of the simulation within the boundries and the increments of the lines to render
#define Nx 100
#define Ny 100
#define lineInc 10

// define the placement and magnitude of the point charge using points as opposed to spacial position
#define chargeX 15
#define chargeY 30
#define chargeMag (1.0e-6)  // 1 microcoulomb

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
 *  @note The function currently assumes there is only one charge present at (0,0) with magnitude e based on the constants defined above
 */
double calcVolt(int Px, int Py) {
    
    double xGridPoint = Xmin + Px * dx;
    double yGridPoint = Ymin + Py * dy;
    
    double distanceToCharge = sqrt(pow(xGridPoint - chargeX, 2) + pow(yGridPoint - chargeY, 2));
    
    //printf("Grid Point (%d, %d) -> (%.2f, %.2f)\n", Px, Py, xGridPoint, yGridPoint); <-- Save for printing in case of error
    
    return (coul * chargeMag) / (distanceToCharge + 1e-9);
}

/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calcVolt for every point in the space
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setupGrid(void) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            V[i][j] = calcVolt(i, j); // Use to calculate the difference
        }
    }
}

/**
 *  @brief calculates either the distance to the wall or the direction of the furthest wall based on the given point
 *
 *  @param Cx The point-x position of the point of reference
 *  @param Cy The point-y postition of the point of reference
 *  @param indexReq The index requested boolean. decides which value is returned
 *
 *  @return Based on indexReq, the function either returns the distance to the furthest wall or the direction of the furthest wall
 *
 *  @note for directions, the index of the array is used as the return value. The following is true:
 *      - 0: Upwards
 *      - 1: Downwards
 *      - 2: Right
 *      - 3: Left
 */
int findFurthestWall(int Cx, int Cy, int indexReq) { // find the wall furthest from the charge
    int distances[4];
    
    // calculate distances toward each wall
    int toRightWall = (Nx - 1 - Cx);
    int toLeftWall = Cx;
    int toTopWall = Cy;
    int toBottomWall = (Ny - 1 - Cy);
    
    // assign distance values to index of direction
    distances[0] = toTopWall;
    distances[1] = toBottomWall;
    distances[2] = toRightWall;
    distances[3] = toLeftWall;
    
    // use linear search to find index of furthest wall
    int maxIndex = 0;
    for (int i = 1; i <= 3; i++) {
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
    
    int lineCount = distanceToWall / lineInc;
    
    // check line count value before calculation and allocation
    if (lineCount <= 0) {
        printf("Error: Not enough space to compute equipotential lines.\n");
        return;
    }
    
    double *lineData = (double *)malloc(lineCount * sizeof(double));
    if (lineData == NULL) {
        fprintf(stderr, "Memory allocation failed for lineData\n");
        return;
    }
    
    for (int i = 0; i < lineCount; i++) {
        switch (directionToFurthestWall) { // use direction to fill in calues every 25 pts
            case 0: // up
                lineData[i] = V[Cx][Cy + (i * lineInc)];
                break;
            case 1: // down
                lineData[i] = V[Cx][Cy - (i * lineInc)];
                break;
            case 2: // right
                lineData[i] = V[Cx + (i * lineInc)][Cy];
                break;
            case 3: // left
                lineData[i] = V[Cx - (i * lineInc)][Cy];
                break;
            default:
                printf("Error analyzing value of direction %d", directionToFurthestWall);
                break;
        }
    }
    
    printf("Progress check.\nLineData:\n");
    for (int i = 0; i < lineCount; i++) {
        printf("\n  %d: %.2e", i, lineData[i]);
    }
    printf("\nCharge (%d, %d) is furthest from wall %d, and is %d points away.\n", chargeX, chargeY, directionToFurthestWall, distanceToWall);
    
    free(lineData);
}

int main(int argc, const char * argv[]) {
    
    setupGrid();
    
    renderEqLines(chargeX, chargeY);
    
    return 0;
}

//
//  main.c
//  Equipotential-Visualizer
//
//  Created by Kenneth Anderson on 2/15/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// define the spacial boundaries of the simulation
#define Xmin -10.0
#define Xmax 10.0
#define Ymin -5.0
#define Ymax 5.0

// define the resolution (points) of the simulation within the boundaries and the increments of the lines to render
#define Nx 2000
#define Ny 1000

// Line increment values
#define lineInc 50 // used in the ray-based rendering (Legacy)
#define voltInc 15 // used in the newer, voltage based rendering

// define coulombs constant
#define coul (8.99e9)

struct charge {
    double spaceX, spaceY, magnitude;
    int pointX, pointY;
};

double V[Nx][Ny];
double dx = (Xmax - Xmin) / (Nx - 1);
double dy = (Ymax - Ymin) / (Ny - 1);

/**
 *  @brief Function calculates the potential difference of any provided point
 *
 *  The function uses the universal formula and currently computes that based on the points distance to the charge x and y
 *
 *  @param charge structure data with the charge values
 *  @param chargeCount total amount of charges
 *  @param xGrid The x grid point of the calculation
 *  @param yGrid The y grid point of the calculation
 *
 *  @return Func returns the value of the potential difference at the given point
 */
double calcVolt(struct charge charge[], int chargeCount, int xGrid, int yGrid) {
    double xPhys = Xmin + xGrid * dx;
    double yPhys = Ymin + yGrid * dy;
    double potential = 0.0;
    for (int i = 0; i < chargeCount; i++) {
        double distanceToCharge = sqrt(pow(xPhys - charge[i].spaceX, 2) + pow(yPhys - charge[i].spaceY, 2));
        potential += (coul * charge[i].magnitude) / (distanceToCharge + 1e-9);
    }
    return potential;
}

/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calcVolt for every point in the space
 *  @param charge the structure containing the data for the first charge
 *  @param chargeCount the total amounts of charges
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setupGrid(struct charge charge[], int chargeCount) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            V[i][j] = calcVolt(charge, chargeCount, i, j); // Use to calculate the difference
            if (i % 500 == 0 && j % 250 == 0) {  // Sample every 500x250 points
                printf("V[%d][%d] = %.2e\n", i, j, V[i][j]);
            }
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
 * @brief function finds the minimum and maximum voltage in the space, and returns either min or max
 *
 * @param type indicates return; 1 for max, 0 for min
 */
double findVoltageExtrema(int type) {
    double Vmin = V[0][0], Vmax = V[0][0]; // create min and max variables
    // use a for loop to calculate values
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            if (V[i][j] < Vmin) Vmin = V[i][j];
            if (V[i][j] > Vmax) Vmax = V[i][j];
            //printf("%.2f is min\n %.2f is max\n\n", Vmin, Vmax);
        }
    }
    return type ? Vmax : Vmin; // return value based on type
}
/**
 * @brief This function approaches the rendering by using the voltage of the space
 *
 * The function performs as follows:
 *  - The function first determines the maximum and minimum voltage in the space
 *  - Then a predetermined voltage incrementation is used to generate equipotential contours
 *  - For every point (i, j), check if V[i][j] is close to any of those threshold voltages, and render each pixel white or black
 *
 *  @param charge the structure containing the data for the first charge
 *  @param chargeCount the total amounts of charges
 */
void renderVoltageBased(struct charge charge[], int chargeCount) {
    printf("Starting beginning extrema sequence...\n");
    
    int numLines = 15;

    // Find global voltage extrema across the field
    double Vmin = findVoltageExtrema(0);
    double Vmax = findVoltageExtrema(1);
    double invMin = 1.0 / Vmax; // corresponds to high voltage (closer to charge)
    double invMax = 1.0 / Vmin; // corresponds to low voltage (farther out)
    
    if (Vmin < 1e-6) Vmin = 1e-6;

    // Avoid log(0) or negative values
    //double logMin = log10(Vmin + 1e-9);
    //double logMax = log10(Vmax);
    
    printf("Sequence succeeded. Rendering %d voltage lines...\n", numLines);

    // Create the image
    unsigned char img[Ny][Nx];
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            img[j][i] = 0;  // Black background
        }
    }

    // Loop over log-spaced voltage targets
    for (int k = 0; k < numLines; k++) {
            double invTarget = invMin + (invMax - invMin) * k / (numLines - 1);
            double target = 1.0 / invTarget;
            double tol = 0.005 * target; // or tune this a bit

            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    double Vhere = V[i][j];
                    if (fabs(Vhere - target) < tol) {
                        img[j][i] = 255;
                    }
                }
            }
        }

    // Set the charge position to a white dot
    for (int i = 0; i < chargeCount; i++) {
        int Cx = charge[i].pointX;
        int Cy = charge[i].pointY;
        img[Cy][Cx] = 255;  // Mark each charge
    }

    // Write the image to disk
    stbi_write_png("./eqLinesRendered.png", Nx, Ny, 1, img, Nx);
}

int main(int argc, const char * argv[]) {
    
    // create array of 4 charges
    struct charge charge[4];
    
    // Create an int to define the amount of charges the user has requested
    int chargeCount;
    
    // create an input variable
    char strIn[2];
    
    printf("Enter the amount of charges: ");
    scanf("%d", &chargeCount);
    
    if(chargeCount < 5 && chargeCount > 0) {
        for(int i = 0; i < chargeCount; i++) {
            // display charge number and request values
            printf("Charge %d:\nEnter the charge X: ", i + 1);
            scanf("%lf", &charge[i].spaceX);
            printf("Enter the charge Y: ");
            scanf("%lf", &charge[i].spaceY);
            printf("Enter the charge size (enter x to default to one positive microcoulomb and -x to default to one negative microcoulomb): ");
            scanf("%s", strIn);
            
            // compare string and assign auto-values as needed
            if(strcmp(strIn, "x") == 0) {
                charge[i].magnitude = (1.0e-6); // one +microcoulomb
            } else if (strcmp(strIn, "-x") == 0) {
                charge[i].magnitude = (-1.0e-6); // 1 -microcoulomb
            } else {
                charge[i].magnitude = atof(strIn); // input
            }
            
            // assign and set the point coordinates of the charges based on the spatial coordinates
            charge[i].pointX = (int)((charge[i].spaceX - Xmin) / dx + 0.5);
            charge[i].pointY = (int)((charge[i].spaceY - Ymin) / dy + 0.5);
        }
    } else {
        printf("Next time, choose a proper value. Terminating sesh.");
        return 0;
    }
    
    // calculate all ∆V values
    printf("Setting up...\n");
    setupGrid(charge, chargeCount);
    printf("Setup Complete.\n\n");
    //setupGrid(charge1.spaceX, charge1.spaceY, charge1.magnitude);
    
    // render the image
    printf("Rendering...\n");
    renderVoltageBased(charge, chargeCount);
    printf("Rendering Complete.\n\n");
    //renderVoltageBased(charge1.pointX, charge1.pointY);
    
    // exit the program
    return 0;
}

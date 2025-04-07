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
#define lineInc 50
#define voltInc 15 // define the voltage increment (15 volts for now I'll tweak as needed)

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
 *  @param Px is the given point x value of the place to compute the difference
 *  @param Py is the given point y value of the place to compute the difference
 *  @param Cx is the charge spacial x value
 *  @param Cy is the charge spacial y value
 *
 *  @return Func returns the value of the potential difference at the given point
 */
double calcVolt(int Px, int Py, double Cx, double Cy, double chargeSize) {
    
    double xGridPoint = Xmin + Px * dx;
    double yGridPoint = Ymin + Py * dy;

    double distanceToCharge = sqrt(pow(xGridPoint - Cx, 2) + pow(yGridPoint - Cy, 2));
    
    //printf("Grid Point (%d, %d) -> (%.2f, %.2f)\n", Px, Py, xGridPoint, yGridPoint); <-- Save for printing in case of error
    
    return (coul * chargeSize) / (distanceToCharge + 1e-9);
}

/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calcVolt for every point in the space
 *  @param Cx the spatial x coordinate of the charge
 *  @param Cy the spatial y coordinate of the charge
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setupGrid(double Cx, double Cy, double chargeSize) {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            V[i][j] = calcVolt(i, j, Cx, Cy, chargeSize); // Use to calculate the difference
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
 * @param Cx the point-x location of the charge
 * @param Cy the point-y position of the charge
 */
void renderVoltageBased(int Cx, int Cy) {
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
            double tol = 0.0005 * target; // or tune this a bit

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
    img[Cy][Cx] = 255;

    // Write the image to disk
    stbi_write_png("./eqLinesInput.png", Nx, Ny, 1, img, Nx);
}

/**
 *  @brief The function creates a png image of the rendered eq lines
 *
 *  The function performs as follows:
 *  - The function first determines which wall (up, down, right, or left) is the furthest from the charge. It also calculates how far away this wall is, in terms of grid points.
 *  - Based on the distance to the furthest wall and a predefined increment (like 10 grid points per step), it calculates how many "rings" or equipotential lines need to be drawn. Each "ring" represents a set of points where the potential difference (∆V) is constant.
 *  - It dynamically allocates memory for an array where each element will represent one of these rings. The size of this array corresponds to the number of rings calculated from the distance to the furthest wall divided by the increment.
 *  - For each ring, the function computes the potential at specific points along the line from the charge to the wall at intervals defined by the increment. These potential values are stored in the array, representing the voltage at each ring's position.
 *  - The function then generates a PNG image where points with a potential that matches one of the calculated ring values (within some tolerance) are colored white. These points form the equipotential lines or "rings". All other points are colored black.
 *
 *  @param Cx the point-x location of the charge
 *  @param Cy the point-y position of the charge
 */
void renderEqLines(int Cx, int Cy) {

    /* NOTE: The "Distance to wall" method will not work when implementing multiple charges. 
    Instead, it may be easier to calculate a potential difference a fixed distance away, 
    and then increment that value for rendering all of the lines instead of calculating the difference of certain points that fall on a ray.
    */

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
    
    double *lineData = (double *)malloc(lineCount * sizeof(double)); // allocate and double check memory allocation
    if (lineData == NULL) {
        fprintf(stderr, "Memory allocation failed for lineData\n");
        return;
    }
    
    // using a for loop and a switch case to calculate the ∆V values to render based on direction to wall
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
            default: // in case the direction function returned something weird (shouldnt happen btw)
                printf("Error analyzing value of direction %d", directionToFurthestWall);
                break;
        }
    }
    
    // The following simply prints some values for double checking accuracy within the terminal
    printf("Progress check.\nLineData:\n");
    for (int i = 0; i < lineCount; i++) {
        printf("\n  %d: %.2e", i, lineData[i]);
    }
    printf("\nCharge (%d, %d) is furthest from wall %d, and is %d points away.\n", Cx, Cy, directionToFurthestWall, distanceToWall);
    
    // create a base for the grace value to get consistent lines
    double graceBase = 250;
    
    // create the img array in B&W
    unsigned char img[Ny][Nx];
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            img[j][i] = 0;  // Black background
        }
    }
    
    // Fill in the img array with color values based on whether they match the ones saved prior
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < lineCount; k++) {
                double grace = graceBase / ((k + 1) * (k + 1));
                if(fabs(V[i][j] - lineData[k]) < grace) {
                    img[j][i] = 255;  // White
                }
            }
        }
    }
    
    // set the charge point to a white pixel
    img[Cy][Cx] = 255; // White
    
    // send the image (Replace with desired file location)
    stbi_write_png("./eqLinesInput.png", Nx, Ny, 1, img, Nx);
    free(lineData);
}

int main(int argc, const char * argv[]) {
    
    struct charge charge1;
    char strIn[2];
    
    printf("Enter the charge X: ");
    scanf("%lf", &charge1.spaceX);
    printf("Enter the charge Y: ");
    scanf("%lf", &charge1.spaceY);
    printf("Enter charge size (enter x to default to one positive microcoulomb and -x to default to one negative microcoulomb): ");
    scanf("%s", strIn);
    
    if(strcmp(strIn, "x") == 0) {
        charge1.magnitude = (1.0e-6);
    } else if (strcmp(strIn, "-x") == 0) {
        charge1.magnitude = (-1.0e-6);
    } else {
        charge1.magnitude = atof(strIn);
    }
    
    charge1.pointX = (int)((charge1.spaceX - Xmin) / dx + 0.5);
    charge1.pointY = (int)((charge1.spaceY - Ymin) / dy + 0.5);
    
    // calculate all ∆V values
    setupGrid(charge1.spaceX, charge1.spaceY, charge1.magnitude);
    
    // render the image
    renderVoltageBased(charge1.pointX, charge1.pointY);
    
    // exit the program
    return 0;
}

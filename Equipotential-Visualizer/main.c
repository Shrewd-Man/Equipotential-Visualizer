// 
//  main.c
//  Equipotential-Visualizer
//
//  Created by Kenneth Anderson on 2/15/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "config.h"
//#include "charge.h"

struct charge {
    double spaceX, spaceY, magnitude;
    int pointX, pointY;
};

double V[NX][NY];
double dx = (XMAX - XMIN) / (NX - 1);
double dy = (YMAX - YMIN) / (NY - 1);

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
    double xPhys = XMIN + xGrid * dx;
    double yPhys = YMIN + yGrid * dy;
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
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            V[i][j] = calcVolt(charge, chargeCount, i, j); // Use to calculate the difference
            if (i % 500 == 0 && j % 250 == 0) {  // Sample every 500x250 points
                printf("V[%d][%d] = %.2e\n", i, j, V[i][j]);
            }
        }
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
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
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
    unsigned char img[NY][NX];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            img[j][i] = 0;  // Black background
        }
    }

    // Loop over log-spaced voltage targets
    for (int k = 0; k < numLines; k++) {
            double invTarget = invMin + (invMax - invMin) * k / (numLines - 1);
            double target = 1.0 / invTarget;
            double tol = 0.005 * target; // or tune this a bit

            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
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
    stbi_write_png("./eqLinesRendered.png", NX, NY, 1, img, NX);
}

int main(int argc, const char * argv[]) {
    
    // create array of 4 charges
    struct charge charge[4];
    
    // Create an int to define the amount of charges the user has requested
    int chargeCount;
    
    // create an input variable
    char strIn[32];
    
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
            charge[i].pointX = (int)((charge[i].spaceX - XMAX) / dx + 0.5);
            charge[i].pointY = (int)((charge[i].spaceY - YMIN) / dy + 0.5);
        }
    } else {
        printf("Next time, choose a proper value. Terminating sesh.");
        return 0;
    }
    
    // calculate all ∆V values
    printf("Setting up...\n");
    setupGrid(charge, chargeCount);
    printf("Setup Complete.\n\n");
    
    // render the image
    printf("Rendering...\n");
    renderVoltageBased(charge, chargeCount);
    printf("Rendering Complete.\n\n");
    
    // exit the program
    return 0;
}

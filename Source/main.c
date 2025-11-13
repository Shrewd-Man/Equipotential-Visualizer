// 
//  main.c
//  Equipotential-Visualizer
//
//  Created by Kenneth Anderson on 2/15/25.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "charge.h"
#include "render.h"
#include "grid.h"
#include "spaceCalc.h"

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
/*double calcVolt(charge_t charge[], int chargeCount, int xGrid, int yGrid) {
    double xPhys = XMIN + xGrid * dx;
    double yPhys = YMIN + yGrid * dy;
    double potential = 0.0;
    for (int i = 0; i < chargeCount; i++) {
        double distanceToCharge = sqrt(pow(xPhys - charge[i].spaceX, 2) + pow(yPhys - charge[i].spaceY, 2));
        if(charge[i].sign == 1){
            potential += (COULOMBS_CONSTANT * charge[i].magnitude) / (distanceToCharge + 1e-9);
        } else {
            potential += (COULOMBS_CONSTANT * (0 - charge[i].magnitude)) / (distanceToCharge + 1e-9);
        }
        //potential += (COULOMBS_CONSTANT * charge[i].magnitude) / (distanceToCharge + 1e-9);
    }
    return potential;
}*/

/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calcVolt for every point in the space
 *  @param charge the structure containing the data for the first charge
 *  @param chargeCount the total amounts of charges
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
/*void setupGrid(charge_t charge[], int chargeCount) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            V[i][j] = calcVolt(charge, chargeCount, i, j); // Use to calculate the difference
            if (i % 500 == 0 && j % 250 == 0) {  // Sample every 500x250 points
                printf("V[%d][%d] = %.2e\n", i, j, V[i][j]);
            }
        }
    }
}*/


int main(int argc, const char * argv[]) {
    
    // create array of 4 charges
    charge_t charge[4];
    
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
            printf("Enter the charge size (enter x to default to one positive microcoulomb and z to default to one negative microcoulomb): ");
            scanf("%s", strIn);
            
            // compare string and assign auto-values as needed
            if(strcmp(strIn, "x") == 0) {
                charge[i].magnitude = (1.0e-6); // 1 microcoulomb
                charge[i].sign = POSITIVE;
            } else if (strcmp(strIn, "z") == 0) {
                charge[i].magnitude = (1.0e-6); // 1 microcoulomb
                charge[i].sign = NEGATIVE;
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
    setupGrid(charge, chargeCount, V);
    printf("Setup Complete.\n\n");
    
    // render the image
    printf("Rendering...\n");
    renderVoltageBased(charge, chargeCount, V);
    printf("Rendering Complete.\n\n");
    
    // exit the program
    return 0;
}

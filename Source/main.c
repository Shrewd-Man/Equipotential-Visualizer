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
            printf("Charge %d:\nEnter the x position of the charge: ", i + 1);
            scanf("%lf", &charge[i].spaceX);
            printf("Enter the y position of charge: ");
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
        printf("Inputted value exceeds maximum. Terminating session.");
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

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
#include "space_calc.h"

/** @brief Voltage at every grid point, indexed V[x][y] */
double V[NX][NY];
/** @brief Grid spacing along x in meters, derived from XMIN, XMAX and NX */
double dx = (XMAX - XMIN) / (NX - 1);
/** @brief Grid spacing along y in meters, derived from YMIN, YMAX and NY */
double dy = (YMAX - YMIN) / (NY - 1);

/**
 * @brief Program entry point: reads the charges, computes the voltage field, and renders it to a PNG
 *
 * The program runs in four steps:
 *  - Asks for the number of charges (1 to 4), then each charge's position and size
 *  - Computes the voltage at every grid point with setup_grid()
 *  - Asks for an output file path (including the '.png' extension)
 *  - Draws the equipotential contours with render_voltage_based() and writes the image
 *
 * @param argc number of command-line arguments (unused)
 * @param argv command-line arguments (unused)
 *
 * @return 0 on completion, also 0 if the requested charge count is outside 1 to 4
 */
int main(int argc, const char * argv[]) {
    
    // create array of 4 charges
    charge_t charge[4];
    
    // Create an int to define the amount of charges the user has requested
    int charge_count;
    
    // create input buffers (charge size entry, and the output image path)
    char str_in[32];
    char output_path[256];

    printf("Enter the amount of charges: ");
    scanf("%d", &charge_count);
    
    if(charge_count < 5 && charge_count > 0) {
        for(int i = 0; i < charge_count; i++) {
            // display charge number and request values
            printf("Charge %d:\nEnter the x position of the charge: ", i + 1);
            scanf("%lf", &charge[i].space_x);
            printf("Enter the y position of charge: ");
            scanf("%lf", &charge[i].space_y);
            printf("Enter the charge size (enter x to default to one positive nanocoulomb and z to default to one negative nanocoulomb): ");
            scanf("%31s", str_in); // width limit leaves room for the null terminator

            // compare string and assign auto-values as needed
            if(strcmp(str_in, "x") == 0) {
                charge[i].magnitude = (1.0e-9); // 1 nanocoulomb
                charge[i].sign = POSITIVE;
            } else if (strcmp(str_in, "z") == 0) {
                charge[i].magnitude = (1.0e-9); // 1 nanocoulomb
                charge[i].sign = NEGATIVE;
            } else {
                // custom value in coulombs: the sign comes from the value, magnitude is stored as absolute
                double value = atof(str_in);
                charge[i].sign = (value < 0) ? NEGATIVE : POSITIVE;
                charge[i].magnitude = fabs(value);
            }
            
            // assign and set the point coordinates of the charges based on the spatial coordinates
            charge[i].point_x = (int)(((charge[i].space_x - XMIN) / dx) + 0.5);
            charge[i].point_y = (int)(((charge[i].space_y - YMIN) / dy) + 0.5);
        }
    } else {
        printf("Inputted value exceeds maximum. Terminating session.");
        return 0;
    }
    
    // calculate all ∆V values
    printf("Setting up...\n");
    setup_grid(charge, charge_count, V);
    printf("Setup Complete.\n\n");

    printf("\nVoltage calculation successful. Please enter the filepath for the rendered image (include filename and '.png' extension): ");
    scanf("%255s", output_path);

    // render the image
    printf("Rendering...\n");
    render_voltage_based(charge, charge_count, V, output_path);
    printf("Rendering Complete.\n\n");
    
    // exit the program
    return 0;
}

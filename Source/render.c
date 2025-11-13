// render.c
#include <stdio.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "render.h"
#include "charge.h"
#include "config.h"

// --- Voltage Extrema Function ---
/**
 * @brief function finds the minimum and maximum voltage in the space, and returns either min or max
 *
 * @param type indicates return; 1 for max, 0 for min
 * @param V array for voltage calcs
 */
double findVoltageExtrema(int type, double V[NX][NY]) {
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

// --- Rendering Function ---

void renderVoltageBased(charge_t charge[], int chargeCount, double V[NX][NY]) {
    printf("Starting beginning extrema sequence...\n");
    
    int numLines = 15;

    // Find global voltage extrema across the field
    double Vmin = findVoltageExtrema(0, V);
    double Vmax = findVoltageExtrema(1, V);
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
    stbi_write_png("/Users/kenny/Desktop/eqLinesRendered.png", NX, NY, 1, img, NX);
}
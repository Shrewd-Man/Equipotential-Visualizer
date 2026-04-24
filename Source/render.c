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
double findVoltageExtrema(charge_t charge[], int chargeCount, int type, double V[NX][NY]) {
    double Vmin = V[0][0], Vmax = V[0][0]; // create min and max variables
    double distanceToCharge;
    int tooClose = 0;
    // use a for loop to calculate values
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            tooClose = 0;

            for (int k = 0; k < chargeCount; k++) {
                distanceToCharge = sqrt(pow((i - charge[k].pointX), 2) + pow((j - charge[k].pointY), 2));
                if (distanceToCharge < 30) {
                    tooClose = 1;
                }
            }
            if (tooClose == 0) {
                if (V[i][j] < Vmin) Vmin = V[i][j];
                if (V[i][j] > Vmax) Vmax = V[i][j];
            }
            //if (V[i][j] < Vmin) Vmin = V[i][j];
            //if (V[i][j] > Vmax) Vmax = V[i][j];
            //printf("%.2f is min\n %.2f is max\n\n", Vmin, Vmax);
        }
    }
    return type ? Vmax : Vmin; // return value based on type
}

void drawContoursForSign(double V[NX][NY], unsigned char img[NY][NX], int sign, double absMax, int numLines) {
    double Vref = (sign > 0) ? fmax(absMax, 1e-12) : fmax(-absMax, 1e-12);
    if (Vref < 1e-12) return;

    double logMin = log10(1e-12);
    double logMax = log10(Vref);

    double logSpacing = (numLines > 1) ? (logMax - logMin) / (numLines - 1) : 0;
    double tol = logSpacing * 0.25;  // Tolerance as fraction of log spacing

    for (int k = 0; k < numLines; k++) {
        double t = (numLines == 1) ? 0.0 : (double)k / (numLines - 1);
        double logTarget = logMin + (logMax - logMin) * t;
        double targetMag = pow(10.0, logTarget);

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double Vhere = V[i][j];
                if ((Vhere > 0 && sign > 0) || (Vhere < 0 && sign < 0)) {
                    double logVhere = log10(fabs(Vhere));
                    if (fabs(logVhere - logTarget) < tol) {
                        img[j][i] = 255;
                    }
                }
            }
        }
    }
}

// --- Rendering Function ---

void renderVoltageBased(charge_t charge[], int chargeCount, double V[NX][NY], char strIn[32]) {
    printf("Starting beginning extrema sequence...\n");
    
    int numLines = 15;

    // Find global voltage extrema across the field
    double Vmin = findVoltageExtrema(charge, chargeCount, 0, V);
    double Vmax = findVoltageExtrema(charge, chargeCount, 1, V);


    double absMin = fmax(fabs(Vmin), fabs(Vmax));
    if (absMin < 1e-12) absMin = 1e-12;

    double logMin = log10(1e-12);  // floor
    double logMax = log10(absMin);


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
    /*for (int k = 0; k < numLines; k++) {
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
    }*/

    if (Vmax > 1e-12) {
        drawContoursForSign(V, img, +1, Vmax, numLines);
    }

    // Draw negative contours
    if (Vmin < -1e-12) {
        drawContoursForSign(V, img, -1, Vmin, numLines);
    }

    // Set the charge position to a white dot
    for (int i = 0; i < chargeCount; i++) {
        int Cx = charge[i].pointX;
        int Cy = charge[i].pointY;
        img[Cy][Cx] = 255;  // Mark each charge
    }

    // Write the image to disk
    stbi_write_png(strIn, NX, NY, 1, img, NX);
}
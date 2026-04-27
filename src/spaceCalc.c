#include <stdio.h>
#include "spaceCalc.h"
#include "charge.h"
#include "config.h"
#include "grid.h"

// --- Single Point Voltage Calculation Algorithm ---
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
double calcVolt(charge_t charge[], int chargeCount, int xGrid, int yGrid) {
    double xPhys = XMIN + xGrid * dx;
    double yPhys = YMIN + yGrid * dy;
    double potential = 0.0;
    for (int i = 0; i < chargeCount; i++) {
        double distanceToCharge = sqrt(pow(xPhys - charge[i].spaceX, 2) + pow(yPhys - charge[i].spaceY, 2));
        double q = (charge[i].sign == POSITIVE) ? charge[i].magnitude : -charge[i].magnitude;
        potential += (COULOMBS_CONSTANT * q) / (distanceToCharge + 1e-9);
        //potential += (COULOMBS_CONSTANT * charge[i].magnitude) / (distanceToCharge + 1e-9);
    }
    return potential;
}

// --- Loop System for Each Point ---
/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calcVolt for every point in the space
 *  @param charge the structure containing the data for the first charge
 *  @param chargeCount the total amounts of charges
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setupGrid(charge_t charge[], int chargeCount, double V[NX][NY]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            V[i][j] = calcVolt(charge, chargeCount, i, j); // Use to calculate the difference
            if (i % 500 == 0 && j % 250 == 0) {  // Sample every 500x250 points
                printf("V[%d][%d] = %.2e\n", i, j, V[i][j]);
            }
        }
    }
}
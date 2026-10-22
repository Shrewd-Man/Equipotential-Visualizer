#include <stdio.h>
#include "space_calc.h"
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
 *  @param charge_count total amount of charges
 *  @param x_grid The x grid point of the calculation
 *  @param y_grid The y grid point of the calculation
 *
 *  @return Func returns the value of the potential difference at the given point
 */
double calc_volt(charge_t charge[], int charge_count, int x_grid, int y_grid) {
    double x_phys = XMIN + x_grid * dx;
    double y_phys = YMIN + y_grid * dy;
    double potential = 0.0;
    for (int i = 0; i < charge_count; i++) {
        double distance_to_charge = sqrt(pow(x_phys - charge[i].space_x, 2) + pow(y_phys - charge[i].space_y, 2));
        double q = (charge[i].sign == POSITIVE) ? charge[i].magnitude : -charge[i].magnitude;
        potential += (COULOMBS_CONSTANT * q) / (distance_to_charge + 1e-9);
        //potential += (COULOMBS_CONSTANT * charge[i].magnitude) / (distance_to_charge + 1e-9);
    }
    return potential;
}

// --- Loop System for Each Point ---
/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calc_volt for every point in the space
 *  @param charge the structure containing the data for the first charge
 *  @param charge_count the total amounts of charges
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setup_grid(charge_t charge[], int charge_count, double V[NX][NY]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            V[i][j] = calc_volt(charge, charge_count, i, j); // Use to calculate the difference
            if (i % 500 == 0 && j % 250 == 0) {  // Sample every 500x250 points
                printf("V[%d][%d] = %.2e\n", i, j, V[i][j]);
            }
        }
    }
}
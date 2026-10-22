// space_calc.h
#ifndef SPACE_CALC_H
#define SPACE_CALC_H

#include "charge.h"
#include "config.h"

// --- VOLTAGE SPACE CALCULATION ALGORITHM ---

/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calc_volt for every point in the space
 *  @param charge the structure containing the data for the first charge
 *  @param charge_count the total amounts of charges
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setup_grid(charge_t charge[], int charge_count, double V[NX][NY]);


#endif
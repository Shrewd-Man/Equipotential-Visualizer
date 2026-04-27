// spaceCalc.h
#ifndef SPACECALC_H
#define SPACECALC_H

#include "charge.h"
#include "config.h"

// --- VOLTAGE SPACE CALCULATION ALGORITHM ---

/**
 *  @brief Function calculates the potential difference of every point in the matrix
 *
 *  The function calls the calcVolt for every point in the space
 *  @param charge the structure containing the data for the first charge
 *  @param chargeCount the total amounts of charges
 *
 *  @note This function is currently a temporary system that does not allow for live rendering
 */
void setupGrid(charge_t charge[], int chargeCount, double V[NX][NY]);


#endif
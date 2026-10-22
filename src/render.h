// render.h
#ifndef RENDER_H
#define RENDER_H

#include "charge.h"
#include "config.h"

// --- VOLTAGE BASED RENDERING ---

/**
 * @brief This function approaches the rendering by using the voltage of the space
 *
 * The function performs as follows:
 *  - The function first determines the maximum and minimum voltage in the space
 *  - Then a predetermined voltage incrementation is used to generate equipotential contours
 *  - For every point (i, j), check if V[i][j] is close to any of those threshold voltages, and render each pixel white or black
 *
 *  @param charge the structure containing the data for the first charge
 *  @param charge_count the total amounts of charges
 *  @param V the voltage coordinate array for which to render
 *  @param output_path the filepath (including the '.png' extension) to write the rendered image to
 */
void render_voltage_based(charge_t charge[], int charge_count, double V[NX][NY], const char *output_path);

#endif
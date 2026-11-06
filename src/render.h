// render.h
#ifndef RENDER_H
#define RENDER_H

#include "charge.h"
#include "config.h"

// --- VOLTAGE BASED RENDERING ---

/**
 * @brief Renders the voltage field as white equipotential contour lines on a black PNG
 *
 * The function performs as follows:
 *  - The function first determines the maximum and minimum voltage in the space
 *  - Then, for the positive and the negative voltages separately, CONTOUR_LINES log-spaced voltage levels are chosen
 *  - A line is drawn (white) wherever two neighbouring points (i, j) fall between different levels, everything else stays black
 *
 *  @param charge the array of charges (each is drawn as a single white dot)
 *  @param charge_count the total amounts of charges
 *  @param V the voltage coordinate array for which to render
 *  @param output_path the filepath (including the '.png' extension) to write the rendered image to
 */
void render_voltage_based(charge_t charge[], int charge_count, double V[NX][NY], const char *output_path);

#endif
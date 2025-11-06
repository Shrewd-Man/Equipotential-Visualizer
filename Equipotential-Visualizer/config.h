// CONFIG.H
#ifndef CONFIG_H
#define CONFIG_H

// --- NECESSARY HEADERS ---
#include <math.h>

// --- SIMULATION PARAMETERS ---

    // Boundaries
#define XMIN -10.0
#define XMAX 10.0
#define YMIN -5.0
#define YMAX 5.0

    // Resolution
#define Nx 2000
#define Ny 1000

    // Voltage-Based Rendering increments
#define voltInc 15

// --- CONSTANTS ---

    // Coulombs Constant
#define coul (8.99e9)


#endif
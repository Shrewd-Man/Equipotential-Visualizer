// CONFIG.H
#ifndef CONFIG_H
#define CONFIG_H

// --- NECESSARY HEADERS ---
#include <math.h>

// --- SIMULATION PARAMETERS ---

    // Boundaries
#define XMIN -10.0 /**< Left edge of the simulated space, in meters */
#define XMAX 10.0  /**< Right edge of the simulated space, in meters */
#define YMIN -5.0  /**< Bottom edge of the simulated space, in meters */
#define YMAX 5.0   /**< Top edge of the simulated space, in meters */

    // Resolution
#define NX 2000 /**< Number of grid points along x (also the output image width in pixels) */
#define NY 1000 /**< Number of grid points along y (also the output image height in pixels) */

    // Voltage-Based Rendering increments
#define VOLT_INC 15 /**< Voltage increment setting (currently unused) */

    // Voltage-Based Rendering contours
#define CONTOUR_LINES 8      /**< Contour lines drawn for each charge sign */
#define CONTOUR_DECADES 2.0  /**< How far below its peak a sign's contours reach, in powers of ten */

// --- CONSTANTS ---

    // Coulombs Constant
#define COULOMBS_CONSTANT (8.99e9) /**< Coulomb's constant k, in N*m^2/C^2 */


#endif
// charge.h

#ifndef CHARGE_H
#define CHARGE_H

// --- Sign Bit ---
/**
 * @brief Polarity of a point charge
 */
typedef enum {
    NEGATIVE = 0, /**< Negative charge */
    POSITIVE = 1  /**< Positive charge */
} sign_t;

// --- Charge Structure ---
/**
 * @brief A single point charge in the simulated space
 *
 * A charge is stored twice: once in physical coordinates (meters) and once as the grid point it
 * lands on. The magnitude is always stored as a positive number; the polarity lives in the sign field.
 */
typedef struct charge_t {
    double space_x;  /**< X position in meters */
    double space_y;  /**< Y position in meters */

    int point_x; /**< X position as a grid point index (0 to NX - 1) */
    int point_y; /**< Y position as a grid point index (0 to NY - 1) */

    double magnitude; /**< Absolute size of the charge in coulombs */

    sign_t sign; /**< Polarity: NEGATIVE (0) or POSITIVE (1) */
} charge_t;

#endif //CHARGE_H
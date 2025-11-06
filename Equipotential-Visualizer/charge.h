// charge.h

#ifndef CHARGE_H
#define CHARGE_H

// --- Sign Bit ---
typedef enum {
    NEGATIVE = 0,
    POSITIVE = 1
} chargeSign;

// --- Charge Structure ---
typedef struct charge {
    double spaceX;  /* Meters */
    double spaceY;  /* Meters */

    int pointX; /* Point-based Coordinates */
    int pointY; /* Point-based Coordinates */

    double magnitude; /* Coulombs */

    chargeSign sign; /* 0 for negative, 1 for positive */
} charge;

#endif CHARGE_H
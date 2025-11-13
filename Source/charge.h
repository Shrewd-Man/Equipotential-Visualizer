// charge.h

#ifndef CHARGE_H
#define CHARGE_H

// --- Sign Bit ---
typedef enum {
    NEGATIVE = 0,
    POSITIVE = 1
} sign_t;

// --- Charge Structure ---
typedef struct charge_t {
    double spaceX;  /* Meters */
    double spaceY;  /* Meters */

    int pointX; /* Point-based Coordinates */
    int pointY; /* Point-based Coordinates */

    double magnitude; /* Coulombs */

    sign_t sign; /* 0 for negative, 1 for positive */
} charge_t;

#endif //CHARGE_H
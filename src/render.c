// render.c
#include <stdio.h>
#include <float.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "render.h"


// --- Voltage Extrema Function ---
/**
 * @brief function finds the minimum and maximum voltage in the space, and returns either min or max
 *
 * @param type indicates return; 1 for max, 0 for min
 * @param V array for voltage calcs
 */
double find_voltage_extrema(charge_t charge[], int charge_count, int type, double V[NX][NY]) {
    // start from the opposite ends so a charge sitting on the V[0][0] corner can't poison the result
    double v_min = DBL_MAX, v_max = -DBL_MAX;
    double distance_to_charge;
    int too_close = 0;
    // use a for loop to calculate values
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            too_close = 0;

            for (int k = 0; k < charge_count; k++) {
                distance_to_charge = sqrt(pow((i - charge[k].point_x), 2) + pow((j - charge[k].point_y), 2));
                if (distance_to_charge < 30) {
                    too_close = 1;
                }
            }
            if (too_close == 0) {
                if (V[i][j] < v_min) v_min = V[i][j];
                if (V[i][j] > v_max) v_max = V[i][j];
            }
            //if (V[i][j] < v_min) v_min = V[i][j];
            //if (V[i][j] > v_max) v_max = V[i][j];
            //printf("%.2f is min\n %.2f is max\n\n", v_min, v_max);
        }
    }
    return type ? v_max : v_min; // return value based on type
}

// --- Contour Drawing Function ---
/**
 * @brief draws the equipotential contour lines for either the positive or the negative potentials
 *
 * The contour levels are spaced logarithmically between the weakest and strongest potential of the
 * requested sign (never more than CONTOUR_DECADES powers of ten below the strongest, so the lines
 * stay in the part of the field that actually has structure).
 *
 * Every point is assigned the number of contour levels it has passed, and a line is drawn wherever
 * two neighbouring points have a different count. This keeps every line a constant, thin width no
 * matter how steep or shallow the potential is, instead of painting a band around each level.
 *
 * @param V the voltage coordinate array
 * @param img the image to draw the lines into (255 for a line)
 * @param sign +1 to draw the contours of the positive potentials, -1 for the negative potentials
 * @param extremum the strongest voltage of this sign (taken outside the charges' exclusion zone)
 * @param num_lines the number of contour lines to draw for this sign (2 or more)
 */
void draw_contours_for_sign(double V[NX][NY], unsigned char img[NY][NX], int sign, double extremum, int num_lines) {
    double v_hi = fabs(extremum);
    if (v_hi < 1e-12 || num_lines < 2 || num_lines > 250) return;

    // find the weakest potential of this sign that exists anywhere in the field
    double v_seen_min = v_hi;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            double mag = sign * V[i][j]; // positive when the point has the requested sign
            if (mag > 0 && mag < v_seen_min) v_seen_min = mag;
        }
    }

    // the dipole's zero line makes the weakest potential tiny, so limit how far down the levels go
    double v_lo = fmax(v_seen_min, fmax(v_hi * pow(10.0, -CONTOUR_DECADES), 1e-12));
    if (v_lo >= v_hi) return;

    // place the levels strictly inside [v_lo, v_hi] so that every one of them is visible
    double log_lo = log10(v_lo);
    double step = (log10(v_hi) - log_lo) / (num_lines + 1);

    // band[i][j] is how many levels the point has passed (0 = below the first level)
    static unsigned char band[NX][NY];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            double mag = sign * V[i][j];
            int b = 0;
            if (mag > 0) {
                b = (int)floor((log10(mag) - log_lo) / step);
                if (b < 0) b = 0;
                if (b > num_lines) b = num_lines;
            }
            band[i][j] = (unsigned char)b;
        }
    }

    // draw a line wherever neighbouring points fall in different bands
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            if (i + 1 < NX && band[i][j] != band[i + 1][j]) {
                img[j][i] = 255;
                img[j][i + 1] = 255;
            }
            if (j + 1 < NY && band[i][j] != band[i][j + 1]) {
                img[j][i] = 255;
                img[j + 1][i] = 255;
            }
        }
    }
}

// --- Rendering Function ---

void render_voltage_based(charge_t charge[], int charge_count, double V[NX][NY], const char *output_path) {
    printf("Starting beginning extrema sequence...\n");

    // Find global voltage extrema across the field (ignores the points right next to the charges)
    double v_min = find_voltage_extrema(charge, charge_count, 0, V);
    double v_max = find_voltage_extrema(charge, charge_count, 1, V);

    printf("Sequence succeeded. Rendering %d voltage lines per charge sign...\n", CONTOUR_LINES);

    // Create the image
    unsigned char img[NY][NX];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            img[j][i] = 0;  // Black background
        }
    }

    // Draw positive contours
    if (v_max > 1e-12) {
        draw_contours_for_sign(V, img, +1, v_max, CONTOUR_LINES);
    }

    // Draw negative contours
    if (v_min < -1e-12) {
        draw_contours_for_sign(V, img, -1, v_min, CONTOUR_LINES);
    }

    // Set the charge position to a white dot
    for (int i = 0; i < charge_count; i++) {
        int cx = charge[i].point_x;
        int cy = charge[i].point_y;
        img[cy][cx] = 255;  // Mark each charge
    }

    // Write the image to disk
    stbi_write_png(output_path, NX, NY, 1, img, NX);
}

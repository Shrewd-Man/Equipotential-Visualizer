// render.c
#include <stdio.h>

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
    double v_min = V[0][0], v_max = V[0][0]; // create min and max variables
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

void draw_contours_for_sign(double V[NX][NY], unsigned char img[NY][NX], int sign, double abs_max, int num_lines) {
    double v_ref = (sign > 0) ? fmax(abs_max, 1e-12) : fmax(-abs_max, 1e-12);
    if (v_ref < 1e-12) return;

    double log_min = log10(1e-12);
    double log_max = log10(v_ref);

    double log_spacing = (num_lines > 1) ? (log_max - log_min) / (num_lines - 1) : 0;
    double tol = log_spacing * 0.25;  // Tolerance as fraction of log spacing

    for (int k = 0; k < num_lines; k++) {
        double t = (num_lines == 1) ? 0.0 : (double)k / (num_lines - 1);
        double log_target = log_min + (log_max - log_min) * t;
        double target_mag = pow(10.0, log_target);

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double v_here = V[i][j];
                if ((v_here > 0 && sign > 0) || (v_here < 0 && sign < 0)) {
                    double log_v_here = log10(fabs(v_here));
                    if (fabs(log_v_here - log_target) < tol) {
                        img[j][i] = 255;
                    }
                }
            }
        }
    }
}

void draw_gray_map() {

    

}

// --- Rendering Function ---

void render_voltage_based(charge_t charge[], int charge_count, double V[NX][NY], const char *output_path) {
    printf("Starting beginning extrema sequence...\n");
    
    int num_lines = 15;

    // Find global voltage extrema across the field
    double v_min = find_voltage_extrema(charge, charge_count, 0, V);
    double v_max = find_voltage_extrema(charge, charge_count, 1, V);


    double abs_min = fmax(fabs(v_min), fabs(v_max));
    if (abs_min < 1e-12) abs_min = 1e-12;

    double log_min = log10(1e-12);  // floor
    double log_max = log10(abs_min);


    // Avoid log(0) or negative values
    //double log_min = log10(v_min + 1e-9);
    //double log_max = log10(v_max);
    
    printf("Sequence succeeded. Rendering %d voltage lines...\n", num_lines);

    // Create the image
    unsigned char img[NY][NX];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            img[j][i] = 0;  // Black background
        }
    }

    // Loop over log-spaced voltage targets
    /*for (int k = 0; k < num_lines; k++) {
        double inv_target = inv_min + (inv_max - inv_min) * k / (num_lines - 1);
        double target = 1.0 / inv_target;
        double tol = 0.005 * target; // or tune this a bit

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                double v_here = V[i][j];
                if (fabs(v_here - target) < tol) {
                    img[j][i] = 255;
                }
            }
        }
    }*/

    if (v_max > 1e-12) {
        draw_contours_for_sign(V, img, +1, v_max, num_lines);
    }

    // Draw negative contours
    if (v_min < -1e-12) {
        draw_contours_for_sign(V, img, -1, v_min, num_lines);
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
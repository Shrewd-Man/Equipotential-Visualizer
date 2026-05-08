// GrayMap.c 

/*---------------------------------------------------------\
| The purpose of this is to create a new rendering system  |
| that maps out the potential difference using a grayscale |
| system where the higher the |∆V|, the darker the shade   |
\---------------------------------------------------------*/

#include <stdio.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "grayMap.h"
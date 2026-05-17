# Equipotential Simulator/Visualizer
<p align="center">
  <img src="images/eqLines2ChargeTest.png" alt="Equipotential Lines Simulation" width="800">
</p>

**A high-resolution 2D electrostatic equipotential line simulator written in C.**

## Features
- User defined point charges
- High resolution rendering and output (2000x1000 by default)
- Contour drawing using LOG based mathematics
- Written in pure C (only external dependency being STB)

## Demonstrated Skills
- **Low-level programming and mathematics**: This program demonstrates low level mathematics and memory management.
- **C Proficiency**: The lack of external dependecies and pure C showcases a proficiency in a highly-tedious and unforgiving programming language.
  - Modular design with header/source separation, structs, enums.
  - Grid-based field computation, distance calculations, avoiding singularities.
- **Image processing**: Integration with `stb_image_write.h` for PNG output to desired location.

<p align="center">
  <img src="images/R2.png" alt="Equipotential Lines Simulation with 3 points" width="800">
</p>

## Installation & Build

```bash
git clone https://github.com/Shrewd-Man/Equipotential-Visualizer.git
cd Equipotential-Visualizer

# Compile (from src directory or adjust paths)
gcc src/*.c -o equipotential -lm -I src
````
Run `./equipotential`

## Usage Example
- 2 charges: +1 nC at (-2, 0), -1 nC at (2, 0) → classic dipole equipotential pattern.
- Supports up to 4 charges currently (easily extensible).

## Roadmap
- **First Change**: addition of a grayscale heatmap-like render to outperform flat contours.
- Modular input to allow implementation in custom applications.
- Real-time data output formatted for custom rendering programs.
- 3D extention.

## License
Use in compliance with the GPL-3.0 license


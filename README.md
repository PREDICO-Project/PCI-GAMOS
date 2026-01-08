# PCI-GAMOS

This repository provides a GAMOS plugin for the simulation of **X-ray Phase-Contrast Imaging (XPCI)** using wave-optical propagation.  
The plugin extends the standard Monte Carlo particle transport in GAMOS/Geant4 by incorporating **wavefront-based modeling**, enabling the simulation of phase-sensitive imaging techniques beyond pure attenuation contrast.

The code is primarily intended for research and development in X-ray imaging, including Talbot-Lau interferometry and propagation-based phase contrast.

The code inherits a custom plug-in developed in [MIMAC](https://github.com/PREDICO-Project/MIMAC). This plug-in is the **GetImage.cc** plug-in which generates an image through a detector. It is used for Snell-only simulations (without Fresnel propagation).


## Features

- Wavefront-based X-ray imaging within GAMOS
- Simulation of:
  - Propagation based Imaging (PBI)
  - Grating based Imaging (GBI)
    - Absorption gratings
    - Phase gratings
    - Talbot-Lau interferometer
- Fresnel propagation of the complex wavefield
- Phase stepping simulations
- Generation of:
  - Intensity images
  - Phase-sensitive signals (e.g. phase gradient)
- Compatible with standard GAMOS workflows and macros


## Requirements

- **GAMOS** (built on Geant4)
- **Geant4** (compatible with your GAMOS version)
- **FFTW3** (for Fresnel propagation)
- C++ compiler compatible with Geant4
- Linux environment (tested on standard Linux distributions)
- python

## Installation

* Download or clone the repository
* Compile the plug-ins of the folder **plug-ings** while GAMOS is setted up. To compile the plug-in type **make** in the terminal while you are inside the **Physics** and **UserActions** subfolders.
* To include the Physics needed to perform PCI write this line in your main simulation file before the initialization.
  ```
  /gamos/physicsList PhysicsListPC
  ```

* Ensure FFTW3 is installed and available.

* **Don't use another Physic List**.
* data/complex_refractive_index folder contains the refraction index for all elements. This folder must not be changed. If you want to modify the folder you must change the path of the folder in the PhaseContrast.cc file, line 211.

## Basic Usage

The plugin is controlled through **GAMOS input macros**, typically using a dedicated input file (e.g. `GetWavefrontTL.in`).

A minimal workflow is:

1. Define the geometry (source, sample, gratings, detector)
2. Enable the wavefront plugin in the macro
3. Configure wavefront and grating parameters
4. Run GAMOS normally

### Example execution

```bash
gamos main_phase_Snell.in
```

We provide some examples to test and learn how to use the developed plug-ins.

* **main_phase_Snell.in** : This macro is used to run a Simulation of a Sphere made of water. This simulation is a Snell-only simulation.
* **main_phase_propagation.in** : This macro is used to run a Simulation of a Sphere made of water. This simulation perform the Wavefront generation and propagation (Fresnel formalism.)

The **geom** folder contains some examples of worlds and elements (GAMOS solids) to put into the defined worlds.
The **spectra** folder contains some spectra to use in the simulations.
The **inputs** folder contains some macros to perform differents tasks like the Image Generation.

We recommend to not modify the folders distribution.

## Key Parameters

Common parameters configured in the wavefront input file include:

* Pixel size of the wavefront grid (X and Y)
* Number of pixels of the wavefront grid (X and Y)
* Propagation distance
* X-ray energy
* Period of absorption or phase gratings
* Relative grating displacement for phase stepping
* Output directory for results

### Example

```bash
/gamos/setParam GetWavefront:NumPixelsX 250
/gamos/setParam GetWavefront:NumPixelsY 250

/gamos/setParam GetWavefront:gridX 15000
/gamos/setParam GetWavefront:gridY 15000
/gamos/setParam GetWavefront:gridSizeX 0.2*um
/gamos/setParam GetWavefront:gridSizeY 0.2*um

/gamos/setParam GetWavefront:zVirtualPlane 404

/gamos/setParam GetWavefront:PropDistance 8.07

/gamos/setParam GetWavefront:debugMode 0

/gamos/setParam DoTalbot:Talbot 1
/gamos/setParam DoTalbot:PeriodG2 1.0*um
/gamos/setParam DoTalbot:displacement 0. 0. 0.
/gamos/setParam DoTalbot:rotationAxis 0. 0. 1.
/gamos/setParam DoTalbot:rotationAngle 0*deg

/gamos/setParam GetWavefront:ResultsFolder output/TL/
/gamos/setParam GetWavefront:OutputFilename test

```

## Output

The plugin generates wave-optical image data for each phase step:

- 2D intensity images
- Separate simulations for:
  - Reference (without object)
  - Object (with sample)

### Output formats

- MetaImage output files (MHD/RAW)

These outputs are typically post-processed into image stacks for phase-contrast reconstruction. We provide a simple script **energyToStack.py**.


## Post-processing

Raw outputs are usually combined and processed using external tools (Python, MATLAB, etc.) to obtain:

- Transmission images
- Differential phase images
- Dark-field images

We recommend to use the work of [XPCIpy](https://github.com/PREDICO-Project/XPCIpy).

### Typical post-processing steps

1. Stack phase-stepping images
2. Apply phase retrieval algorithms
3. Normalize reference and object data
4. Visualize phase-contrast signals

## Intended Use

This plugin is designed for:

* Method development in X-ray phase-contrast imaging
* Verification and validation of phase-contrast algorithms
* Performance and feasibility studies
* Research applications only (not intended for clinical use)

## Limitations

- Primarily monochromatic simulations
- Assumes coherent or partially coherent illumination depending on configuration

## Citation

If you use this plugin in your work, please cite the corresponding publication describing the method and implementation.

```
@article{SANCHEZLARA2026105716,
title = {X-ray phase contrast imaging in GAMOS},
journal = {Physica Medica},
volume = {142},
pages = {105716},
year = {2026},
issn = {1120-1797},
doi = {https://doi.org/10.1016/j.ejmp.2026.105716},
url = {https://www.sciencedirect.com/science/article/pii/S1120179726000037},
author = {V. Sanchez-Lara and F.R. Lozano and C. Huerga and Luis C. Martinez-Gomez and D. Garcia-Pinto},
keywords = {X-ray, Phase contrast, Wavefront, Monte Carlo, Simulation, Geant4},
}
```

## Contact

For questions, issues or contributions, please open an issue on GitHub or contact to vicsan05@ucm.es.

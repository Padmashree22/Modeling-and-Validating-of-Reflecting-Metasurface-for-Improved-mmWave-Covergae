===========================================
         Ansys HFSS Simulation Project
===========================================

# Topic: Modeling and Validation of Reflecting Metasurfaces for Improved mmWave Coverage
# Master Thesis by Padmashree Reddy Malur Vemana
# Supervisor: Simon Häger
# Professor: Dr.-Ing.Christian Wietfeld

# Faculty of Electrical ENgineering and Information Technology
# Communication Netwoeks Institute
# Technische universität Dortmund

===========================================

## Overview

This repository contains the Ansys HFSS simulation models actively used for the thesis. The `Ansys` folder within the `Codebase` directory includes four main folders:

1. **HELIOS_array**
2. **HELIOS_footprint**
3. **Triangulations_Ansys_models**
4. **Urban_scenario_model**

### HELIOS_array and HELIOS_footprint

These folders contain the main HELIOS model required for the thesis. `HELIOS_footprint` represents the 3D component export of the main model, suitable for creating multiple arrays of one's choice.

### Triangulations_Ansys_models

Contains different models for the case study mentioned in the thesis.

### Urban_scenario_model

Replicates reflectors from the urban scenario mentioned in the thesis.

===========================================

## Importing HELIOS_footprint as 3D Component

1. Open a new Ansys project.
2. Right-click on `HFSSDesign1` and select Solution Type as SBR+.
3. Create a Coordinate System using "offset origin" in the Draw panel.
4. In the left dialogue box extension, find "3D Components," right-click, and select "Browse 3D Component."
5. Choose the location of the 3D component.
6. A pop-up box will appear for variable changes (slope angles: alpha and beta, reflector size: length and width). Enter the correct details.
7. Click enter, and the new model will be imported into your current project.

===========================================

## Running Ansys Simulation without Python Scripts

1. Open the required .aedt file.
2. Under the project name, find `HFSSDesign1(SBR+)` on the left dialogue box, which contains all parameters and variables.
3. Set up slope angles, reflector dimensions, incident angles, and frequency under various tabs.
4. Export RCS Total Far-field reports as 3D Polar Plot in dB under "Results."
5. Set the "Radiation" type to infinite sphere for far-field setup with angular resolution.
6. Save the current work and check for used metrics.
7. Under the "Simulation" panel, validate the model. If no errors, then "Analyze All."

===========================================

## Running Ansys Simulation with Python Scripts

1. Follow the procedure in the second point.
2. Instead of "Validate" and "Analyze All," edit the Python script with the required parameters, folder paths, and Ansys project path.
3. Under the "Tools" panel, select "Run Script" and choose the .py file on the current system.
4. The automation of the simulation will start.


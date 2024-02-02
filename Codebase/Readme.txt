===========================================
             Codebase README
===========================================

# Topic: Modeling and Validation of Reflecting Metasurfaces for Improved mmWave Coverage
# Master Thesis by Padmashree Reddy Malur Vemana
# Supervisor: Simon Häger
# Professor: Dr.-Ing.Christian Wietfeld

# Faculty of Electrical Engineering and Information Technology
# Communication Networks Institute
# Technische universität Dortmund

===========================================
              Codebase Overview
===========================================

The `Codebase` folder is structured into four subfolders, each serving a specific purpose in the thesis project:

1. **Ansys**
   - Location: `Codebase/Ansys/`
   - Description: This folder contains productively working simulation models built using the 4-core Ansys HFSS simulator (.aedt file format). Key components include:
     - **HELIOS Model**: Located in the "HELIOS_array" folder, this model is central to the thesis. The "HELIOS_Footprint" folder contains its footprint.
     - **Triangulation Models**: Found in the "Triangulation_Ansys_models" folder, these models support the triangulation of squared models.
     - **Urban Scenario Model**: Present in the "Urban_scenario_model" folder, this model replicates an urban scenario discussed in the thesis.

2. **csvFiles**
   - Location: `Codebase/csvFiles/`
   - Description: This folder contains various CSV files used for generating meaningful results throughout the thesis. Subfolders include:
     - **HELIO_array_Simulation_csvFiles**: Validating the HELIOS model with array configurations from 1x1 to 8x8.
     - **UrbanScenario_Narrow_Broad_Sim_Analytical**: Simulated and analytically calculated results for an urban scenario with narrow and broad beams.
     - **UrbanScenario_Analytical_Optimized_search_csvFiles**: Exhaustive search data in Matlab to find optimized geometry for the analytical HELIOS model.
     - **Computing_yime_csvFiles**: Contains computing time data for 1-core analytical Matlab, 4-core, and 128-core Ansys simulations.
     - **Casestudy_UrbanScenario_orthogonal_incident_size_vary**: Files for a case study in an urban scenario with orthogonal incident waves and varying reflecting element sizes.
     - **Casestudy_Triangulation_Sim_csvFiles**: Simulation files for square and triangle models in different orientations.
     - **Casestudy_Flatplate_Comparison_csvFiles**: Ansys simulation files for a case study comparing different array flatplates.

3. **matlab**
   - Location: `Codebase/matlab/`
   - Description: This folder contains Matlab codes (.m files) organized into subfolders:
     - **HELIOS Model**: Main analytical modeling code and supporting functions.
     - **IRS Models**: Matlab code for three IRS models presented in the thesis and a comparison script.
     - **Urban Scenario**: Matlab codes implementing HELIOS and IRS models in urban scenarios.
     - **Channel Models**: Implementation of three different channel models and their radar realization.
     - **Case Studies**: Matlab code for various case studies conducted in the thesis.
     - **Computing Time**: Matlab code comparing computing time between 1-core analytical model, 4-core, and 128-core simulation models.

4. **python_scripts**
   - Location: `Codebase/python_scripts/`
   - Description: This folder contains Python scripts (.py files) used by the Ansys simulation software for automation:
     - **PythonScript_code_automation.py**: General code for automating simulations with customizable parameters.
     - **PythonScript_ComputingTimeCalculation.py**: Script for automating simulation tool with multiple parameter changes and computing time storage.

===========================================

===========================================

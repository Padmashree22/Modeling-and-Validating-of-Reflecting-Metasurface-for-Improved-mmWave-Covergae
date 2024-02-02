===========================================
             MATLAB README
===========================================

# Topic: Modeling and Validation of Reflecting Metasurfaces for Improved mmWave Coverage
# Master Thesis by Padmashree Reddy Malur Vemana
# Supervisor: Simon Häger
# Professor: Dr.-Ing.Christian Wietfeld

# Faculty of Electrical Engineering and Information Technology
# Communication Networks Institute
# Technische universität Dortmund

===========================================
              Overview
===========================================

This repository contains a folder named `Codebase`, which encompasses MATLAB codes organized into different topics. The `matlab` folder is further subdivided into six main folders:

1. **Case Studies**
   - Location: `Codebase/matlab/Case Studies/`
   - Description: This folder contains MATLAB codes related to various case studies, each with its own subfolder. Additionally, dependency functions are included within each subfolder.

   - Dependencies:
      - **Flatplate_Comparison**
         - Folder: `Codebase/matlab/Case Studies/Flatplate_Comparison`
         - Dependencies: `csvFiles/Casestudy_Flatplate_Comparison_csvfiles`

      - **Triangulation_study**
         - Folder: `Codebase/matlab/Case Studies/Triangulation_study`
         - Dependencies: `csvFiles/Casestudy_Triangulation_Sim_csvfiles`

      - **Triangulation_study_comparison_square**
         - Folder: `Codebase/matlab/Case Studies/Triangulation_study_comparison_square`
         - Dependencies: `csvFiles/Casestudy_Triangulation_Sim_csvfiles`

      - **UrbanScenario_Casestudy_reflectorsize_vary**
         - Folder: `Codebase/matlab/Case Studies/UrbanScenario_Casestudy_reflectorsize_vary`
         - Dependencies: `csvFiles/Casestudy_UrbanScenario_reflectorsize_vary`

      - **Urbanscenario_size_beamwidth_orthogonal_vary**
         - Folder: `Codebase/matlab/Case Studies/Urbanscenario_size_beamwidth_orthogonal_vary`
         - Dependencies: `csvFiles/Casestudy_urbanScenario_orthogonal_incident_size_vary`

2. **Channel Models**
   - Location: `Codebase/matlab/Channel Models`
   - Description: This folder contains MATLAB codes related to channel models. No additional CSV files are required for these codes.

3. **Computing Time**
   - Location: `Codebase/matlab/Computing Time`
   - Description: MATLAB codes related to computing time analysis.
   - Dependencies: `csvFiles/Computingtime_csvfiles`

4. **HELIOS Model**
   - Location: `Codebase/matlab/HELIOS Model`
   - Description: MATLAB codes related to the HELIOS model.
   - Dependencies: `csvFiles/HELIOS_array_Simulation_csvFiles`

5. **IRS Model**
   - Location: `Codebase/matlab/IRS Model`
   - Description: MATLAB codes related to the IRS model. No additional CSV files are required for these codes.

6. **Urban Scenario**
   - Location: `Codebase/matlab/Urban Scenario`
   - Description: MATLAB codes related to urban scenarios. Each subfolder has its dependency functions.

   - Dependencies:
      - **UrbanScenario_Sim_RCS_compare**
         - Folder: `Codebase/matlab/Urban Scenario/UrbanScenario_Sim_RCS_compare`
         - Dependencies: `csvFiles/UrbanScenario_narrow_broad_Sim_Analytical`

      - **Heatmap_search_alphabeta_with_csvFile**
         - Folder: `Codebase/matlab/Urban Scenario/Heatmap_search_alphabeta_with_csvFile`
         - Dependencies: `csvFiles/UrbanScenario_AnalyticalOptimized_search_csvfiles`

===========================================
               Usage
===========================================

Please download the entire `Codebase/...` folder to run any specific MATLAB code. Make sure to download the respective CSV files mentioned under each topic to execute the MATLAB codes successfully.

===========================================

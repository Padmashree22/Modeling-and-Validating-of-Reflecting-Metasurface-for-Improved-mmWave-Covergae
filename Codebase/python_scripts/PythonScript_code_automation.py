import ScriptEnv
import os
import time

# Path to the current folder [on the simulation computer!]
path_to_current_folder = "C:/Users/path_to_current_folder"
# Name of this HFSS Project
project_name = "project_name"

# Create a folder for simulation results if needed
results_folder_name = project_name + "file_name"
if not os.path.exists(results_folder_name):
    os.makedirs(results_folder_name)

# Interaction Beginning with Ansys
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.SetActiveProject("Model_triangle")
oDesign = oProject.SetActiveDesign("HFSSDesign1")

# Parameters to change
el_angle= [90]  # Add theta_in values
az_angle= [0,10,20,30,40,50,60,70,80,90]  # Add phi_in values

L_phi_in = len(az_angle)  # Length of the phi_in array
L_theta_in = len(el_angle)  # Length of the theta_in array

# Loop through all combinations of variables
for phi_in_index in range(L_phi_in):
    for theta_in_index in range(L_theta_in):
        # Change azimuth, elevation, phi_in, and theta_in angles
        oDesign.ChangeProperty(
            [
                "NAME:AllTabs",
                [
                    "NAME:LocalVariableTab",
                    [
                        "NAME:PropServers",
                        "LocalVariables"
                    ],
                    [
                        "NAME:ChangedProps",
                        [
                            "NAME:az_angle",
                            "Value:=", str(az_angle[phi_in_index]) + "deg"
                        ],
                        [
                            "NAME:el_angle",
                            "Value:=", str(el_angle[theta_in_index]) + "deg"
                        ]
                    ]
                ]
            ]
        )

# Click "Fit All" Button
oDesktop.RestoreWindow()

# Solve via SBR+ Solver
oDesign.AnalyzeAll()

# Save Ansys Project File
oProject.Save()

# Update FarField Pattern Report using Simulation Results
oModule = oDesign.GetModule("ReportSetup")
oModule.UpdateReports(["FarField_Report"])

# Click "Fit All" Button
oDesktop.RestoreWindow()

# Export Reflection Pattern
file_name = "_phi_in=" + str(az_angle[phi_in_index]) + "_theta_in=" + str(el_angle[theta_in_index]) + ".csv"
oModule.ExportToFile("FarField_Report", os.path.join(path_to_current_folder, results_folder_name, file_name))

# Save Ansys Project File
oProject.Save()

# DONE
# Click "Fit All" Button
oDesktop.RestoreWindow()

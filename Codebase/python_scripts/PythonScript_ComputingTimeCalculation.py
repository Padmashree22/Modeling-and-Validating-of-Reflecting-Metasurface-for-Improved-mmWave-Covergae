import ScriptEnv
import os
import time
import csv

# Path to the current folder on the simulation computer
path_to_current_folder = "C:/Users/vemana/Desktop/HELIOS_1109/"
project_name = "2x2_HELIOS"

# Initialize Ansys HFSS
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()

oProject = oDesktop.SetActiveProject("2x2_HELIOS") 
oDesign = oProject.SetActiveDesign("HFSSDesign1")
#oDesign = oProject.GetActiveDesign()


# Define the list of values for each variable
alpha11 = [2,4,6,8]
beta11 = [2,4,6,8]  # Example values for beta11
az_angle = [0,5,10,15,20]  # Example values for azimuth angle
el_angle = [90,95,100,105,110]  # Example values for elevation angle

# Function to run a simulation for a specific set of parameters
def run_simulation(alpha, beta, az, el, results_folder):
    oDesign.ChangeProperty([
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
                    "NAME:alpha11",
                    "Value:=", str(alpha) + "deg"
                ],
                [
                    "NAME:beta11",
                    "Value:=", str(beta) + "deg"
                ],
                [
                    "NAME:az_angle",
                    "Value:=", str(az) + "deg"
                ],
                [
                    "NAME:el_angle",
                    "Value:=", str(el) + "deg"
                ]
            ]
        ]
    ])

	# Perform the simulation and save results
    oDesktop.RestoreWindow()
    # Measure start time
    start_time = time.time()
    oDesign.AnalyzeAll()
    # Measure end time
    end_time = time.time()

    oProject.Save()
    oModule = oDesign.GetModule("ReportSetup")
    oModule.UpdateReports(["FarField_Report"])
    oDesktop.RestoreWindow()
    #file_name = "script1.csv"
    #oModule.ExportToFile("FarField_Report", os.path.join(results_folder, file_name))
    oProject.Save()

    # Calculate the elapsed time
    elapsed_time = end_time - start_time

    return elapsed_time
	# Create a CSV file to store simulation times #(CHANGE 3)
	results_folder_name = project_name + "scripts"
	if not os.path.exists(results_folder_name):
    		os.makedirs(results_folder_name)

csv_file_path = os.path.join(results_folder_name, "128core_HELIOS_2x2_2deg_res.csv")
# Check if the file already exists
file_exists = os.path.exists(csv_file_path)

with open(csv_file_path, mode="a") as csv_file:
    fieldnames = ["alpha11", "beta11", "az_angle", "el_angle", "simulation_time"]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    # Write the header row if the file doesn't exist
    if not file_exists:
        writer.writeheader()

    for alpha in alpha11:
        for beta in beta11:
            for az in az_angle:
                for el in el_angle:
                    elapsed_time = run_simulation(alpha, beta, az, el, results_folder_name)

                    # Write simulation times to the CSV file
                    writer.writerow({
                        "alpha11": alpha,
                        "beta11": beta,
                        "az_angle": az,
                        "el_angle": el,
                        "simulation_time": elapsed_time
                    })

# Close Ansys HFSS
oDesktop.RestoreWindow()
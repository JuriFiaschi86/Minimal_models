#! /usr/bin/python3

# This script runs the micrOMEGA parameter scan
# Several SPheno spectrum files are generated
# The micrOMEGA code is called for each generated spectrum and the results are stored in an output file


import shutil
import subprocess
import os
import stat
import glob
import re
import itertools

import Casas_Ibarra_function

###########################
### 0. PRELIMINARY SETTINGS
###########################
# Setting model names
model_name = "T13Balpha02g"
# Setting codes versions
SPheno_version = "4.0.3"
micrOMEGA_version = "4.3.5"
# Setting main directories
SPHENO_folder = "/home/users/fiaschi/Dokumente/SPheno/"
mircOMEGA_folder = "/home/users/fiaschi/Dokumente/micrOMEGA/"
local_path = os.path.dirname(os.path.abspath(__file__)) + "/"
# Setting son directory
SPHENO_local_folder = SPHENO_folder + "SPheno-" + SPheno_version + "/"
mircOMEGA_local_folder = mircOMEGA_folder + "micromegas_" + micrOMEGA_version + "/"
SPHENO_work_folder = SPHENO_local_folder + model_name + "/"
mircOMEGA_work_folder = mircOMEGA_local_folder + model_name + "/"

# Setting input/output file names
SPheno_input_template_name = "LesHouches.in." + model_name + "_low"
SPheno_new_template_name = "Fixed_LesHouches.in." + model_name
SPheno_input_list = []
SPheno_output_list = []

# Flag for Casas-Ibarra (0 = off / 1 = on)
Casas_Ibarra_flag = 1

# Flag for Simon output format of micrOMEGA (0 = my output format / 1 = Simon output format)
Simon_micrOMEGA = 1


############################
### 1. IDENTIFY MODEL INPUTS
############################
# Identify the model parameters
model_parameters = ["lambdaInput", "mCPsiInput", "mpsipsipInput", "lambda4Input", "lambda5Input"]
model_blocks = ["MPHI2IN", "LAM1IN", "LAM6IN"]

# Setting the fixed parameters
parameter_names_fixed = ["lambdaInput","mCPsiInput","mpsipsipInput","lambda4Input","lambda5Input"]
parameter_table_fixed = [0.50, 1.0E+03, 1.5E+03, 5.0E-03, 5.0E-03]

block_names_fixed = ["MPHI2IN", "LAM1IN", "LAM6IN"]
block_entries_fixed = [(((1.4E+03)**2, 0), (0, (3.0E+03)**2)), ((0, 0), (0, 0)), ((0.0003, 0.0045), (0.0003, 0.0045), (0.003, 0.0045))]

# Setting the parameters of the scan
#parameter_names_scan = ["mCPsiInput","lambda4Input"]
#parameter_table_scan = [(5000, 5500), (0.005, 0.010)]
parameter_names_scan = ["lambda4Input"]
parameter_table_scan = [(1.0E-03, 5.0E-03, 1.0E-02, 5.0E-02, 1.0E-01)]
#parameter_names_scan = []
#parameter_table_scan = []

### Scan over blocks
#block_names_scan = ["MPHI2IN", "LAM1IN"]
#block_entries_scan = [(((30000, 0), (0, 30000)), ((40000, 0), (0, 40000))), (((0.01, 0.01), (0.01, 0.01)), ((0.02, 0.02), (0.02, 0.02)))]
block_names_scan = []
block_entries_scan = []

# Number of generations of scalars (> 1 necessary for neutrino masses)
n_generation = 2

# Name of the necessary parameters for Casas-Ibarra parametrisation
Casas_Ibarra_input = ["mCPsiInput", "mpsipsipInput", "MPHI2IN", "LAM1IN", "lambda4Input", "lambda5Input"]
Casas_Ibarra_output = "LAM6IN"
Casas_Ibarra_input_values = []


##############################################
### 2. DEFINE FUNCTIONS TO WRITE SPHENO INPUTS
##############################################

# Function to write a new input value for a parameter in the Spheno input file
def write_spheno_input_param(input_text, parameter_name, new_value):
    old_parameter_line = re.findall("^.*" + parameter_name + ".*$", input_text, re.MULTILINE)[0]
    if int(old_parameter_line[:4]) < 10:
        space_front = 5
    else:
        space_front = 6
    space_back = 13
    if old_parameter_line[space_front] == "-":
        space_back = space_back + 1
    parameter_value_string_old = old_parameter_line[space_front:][:space_back]
    parameter_value_string_new = "{:.7E}".format(new_value)
    new_parameter_line = old_parameter_line.replace(parameter_value_string_old, parameter_value_string_new)
    input_text = input_text.replace(old_parameter_line, new_parameter_line)
    return input_text


# Function to write a new input values for a block in the Spheno input file
def write_spheno_input_block(input_text, block_name, new_values):
    for i in range(0, len(new_values)):
        for j in range(0, len(new_values[i])):
            parameter_name = block_name.lower()[:-2] + "\(" + str(i+1) + "," + str(j+1) + "\)"
            old_parameter_line = re.findall("^.*" + parameter_name + ".*$", input_text, re.MULTILINE)[0]
            if old_parameter_line[7] == "-":
                space = 14
            else:
                space = 13
            parameter_value_string_old = old_parameter_line[7:][:space]
            parameter_value_string_new = "{:.7E}".format(new_values[i][j])
            new_parameter_line = old_parameter_line.replace(parameter_value_string_old, parameter_value_string_new)
            input_text = input_text.replace(old_parameter_line, new_parameter_line)
    return input_text


# Function to change the flags in the Spheno input file for micrOMEGA run
def write_spheno_flag(input_text):
    old_line_flag_13 = " 13 1               # 3-Body decays: none (0), fermion (1), scalar (2), both (3) "
    new_line_flag_13 = " 13 0               # 3-Body decays: none (0), fermion (1), scalar (2), both (3) "
    old_line_flag_16 = " 16 1              # One-loop decays "
    new_line_flag_16 = " 16 0              # One-loop decays "
    old_line_flag_50 = " 50 1               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) "
    new_line_flag_50 = " 50 0               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) "
    old_line_flag_55 = " 55 0               # Calculate loop corrected masses "
    new_line_flag_55 = " 55 1               # Calculate loop corrected masses "
    old_line_flag_77 = " 77 0               # Output for MicrOmegas (running masses for light quarks; real mixing matrices)   "
    new_line_flag_77 = " 77 1               # Output for MicrOmegas (running masses for light quarks; real mixing matrices)   "
    input_text = input_text.replace(old_line_flag_13, new_line_flag_13)
    input_text = input_text.replace(old_line_flag_16, new_line_flag_16)
    input_text = input_text.replace(old_line_flag_50, new_line_flag_50)
    input_text = input_text.replace(old_line_flag_55, new_line_flag_55)
    input_text = input_text.replace(old_line_flag_77, new_line_flag_77)
    return input_text

# Function to obtain the Casas-Ibarra output block LAM6IN
def write_Casas_Ibarra_parametrisation(input_text):
    for i in range(0, len(Casas_Ibarra_input)):
        if Casas_Ibarra_input[i] == "MPHI2IN" or Casas_Ibarra_input[i] == "LAM1IN":
            parameter_block = []
            for n in range(0, n_generation):
                for m in range(0, n_generation):
                    parameter_name = Casas_Ibarra_input[i].lower()[:-2] + "\(" + str(n+1) + "," + str(m+1) + "\)"
                    old_parameter_line = re.findall("^.*" + parameter_name + ".*$", input_text, re.MULTILINE)[0]
                    parameter_block.append(float(old_parameter_line[7:][:13]))
                    
            parameter_matrix = ((parameter_block[0], parameter_block[1]), (parameter_block[2], parameter_block[3]))
            Casas_Ibarra_input_values.append(parameter_matrix)
            
        else:
            parameter_name = Casas_Ibarra_input[i]
            old_parameter_line = re.findall("^.*" + parameter_name + ".*$", input_text, re.MULTILINE)[0]
            if int(old_parameter_line[:4]) < 10:
                space = 5
            else:
                space = 6
            Casas_Ibarra_input_values.append(float(old_parameter_line[space:][:13]))
    
    output = Casas_Ibarra_function.Casas_Ibarra(Casas_Ibarra_input_values)
    Casas_Ibarra_output_value = ((output[0,0], output[0,1]), (output[1,0], output[1,1]), (output[2,0], output[2,1]))
    return Casas_Ibarra_output_value


###############################
### 3. WRITE SPHENO INPUT FILES
###############################

# Open template Spheno input file
SPheno_input_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_input_template_name, "r")
text = SPheno_input_template.read()
SPheno_input_template.close()

# Write the flags
text = write_spheno_flag(text)
# Write fixed parameters
for k in range(0, len(parameter_names_fixed)):
    text = write_spheno_input_param(text, parameter_names_fixed[k], parameter_table_fixed[k])
# Write fixed blocks
for m in range(0, len(block_names_fixed)):
    text = write_spheno_input_block(text, block_names_fixed[m], block_entries_fixed[m])

### Create a Spheno input file for each combination of the parameters to scan

if len(parameter_names_scan) and len(block_names_scan):
# Case of mixed scan on parameters and blocks
    combine_param = list(itertools.product(*[list(itertools.product(*parameter_table_scan)), list(itertools.product(*block_entries_scan))]))
    
    for i in range(0, len(combine_param)):
        SPheno_input_new_name = ""
        for n in range(0, len(combine_param[i][0])):
            text = write_spheno_input_param(text, parameter_names_scan[n], combine_param[i][0][n])
            SPheno_input_new_name = SPheno_input_new_name + parameter_names_scan[n] + "_" + str(combine_param[i][0][n]) + "_" 
        for m in range(0, len(combine_param[i][1])):
            text = write_spheno_input_block(text, block_names_scan[m], combine_param[i][1][m])
            SPheno_input_new_name = SPheno_input_new_name + block_names_scan[m] + "_" + str(combine_param[i][1][m][0][0]) + "_" 
        
        if Casas_Ibarra_flag:
            Casas_Ibarra_output_value = write_Casas_Ibarra_parametrisation(text)
            text = write_spheno_input_block(text, Casas_Ibarra_output, Casas_Ibarra_output_value)
        
        SPheno_input_new_name = SPheno_input_new_name[:-1] + ".in"
        SPheno_input_new = open(SPHENO_folder + SPheno_input_new_name, "w")
        SPheno_input_new.write(text)
        SPheno_input_new.close()
        SPheno_input_list.append(SPheno_input_new_name)


elif len(parameter_names_scan):
# Case of scan on parameters only
    combine_param = list(itertools.product(*parameter_table_scan))
    
    for i in range(0, len(combine_param)):
        SPheno_input_new_name = ""
        for j in range(0, len(parameter_names_scan)):
            text = write_spheno_input_param(text, parameter_names_scan[j], combine_param[i][j])
            SPheno_input_new_name = SPheno_input_new_name + parameter_names_scan[j] + "_" + str(combine_param[i][j]) + "_" 
        
        if Casas_Ibarra_flag:
            Casas_Ibarra_output_value = write_Casas_Ibarra_parametrisation(text)
            text = write_spheno_input_block(text, Casas_Ibarra_output, Casas_Ibarra_output_value)
        
        SPheno_input_new_name = SPheno_input_new_name[:-1] + ".in"
        SPheno_input_new = open(SPHENO_folder + SPheno_input_new_name, "w")
        SPheno_input_new.write(text)
        SPheno_input_new.close()
        SPheno_input_list.append(SPheno_input_new_name)

        
elif len(block_names_scan):
# Case of scan on blocks only
    combine_param = list(itertools.product(*block_entries_scan))
    
    for i in range(0, len(combine_param)):
        SPheno_input_new_name = ""
        
        for m in range(0, len(combine_param[i])):
            text = write_spheno_input_block(text, block_names_scan[m], combine_param[i][m])
            SPheno_input_new_name = SPheno_input_new_name + block_names_scan[m] + "_" + str(combine_param[i][m][0][0]) + "_" 
        
        if Casas_Ibarra_flag:
            Casas_Ibarra_output_value = write_Casas_Ibarra_parametrisation(text)
            text = write_spheno_input_block(text, Casas_Ibarra_output, Casas_Ibarra_output_value)
        
        SPheno_input_new_name = SPheno_input_new_name[:-1] + ".in"
        SPheno_input_new = open(SPHENO_folder + SPheno_input_new_name, "w")
        SPheno_input_new.write(text)
        SPheno_input_new.close()
        SPheno_input_list.append(SPheno_input_new_name)

else:
# Case of single parameter point (no scan)
    SPheno_input_new_name = ""
    #for i in range(0, len(parameter_names_fixed)):
        #SPheno_input_new_name = SPheno_input_new_name + parameter_names_fixed[i] + "_" + str(parameter_table_fixed[i]) + "_" 
    #for i in range(0, len(block_names_fixed)-1):
        #SPheno_input_new_name = SPheno_input_new_name + block_names_fixed[i] + "_" + str(block_entries_fixed[i][0][0]) + "_" 
    
    SPheno_input_new_name = model_name + "_"
    
    if Casas_Ibarra_flag:
        Casas_Ibarra_output_value = write_Casas_Ibarra_parametrisation(text)
        text = write_spheno_input_block(text, Casas_Ibarra_output, Casas_Ibarra_output_value)
    
    SPheno_input_new_name = SPheno_input_new_name[:-1] + ".in"
    SPheno_input_new = open(SPHENO_folder + SPheno_input_new_name, "w")
    SPheno_input_new.write(text)
    SPheno_input_new.close()
    SPheno_input_list.append(SPheno_input_new_name)


#################
### 3. RUN SPHENO
#################
# Run SPheno for each input files
for SPheno_input_new_name in SPheno_input_list:
    print("\n")
    print("Producing SPheno spectrum for " + SPheno_input_new_name + "\n")
    SPheno_output_new_name = SPheno_input_new_name[:-3] + ".out"
    SPheno_output_list.append(SPheno_output_new_name)    
    os.chdir(SPHENO_folder)
    command = SPHENO_folder + "SPheno-" + SPheno_version +"/bin/SPheno" + model_name + " " + SPHENO_folder + SPheno_input_new_name + " " + SPHENO_folder + SPheno_output_new_name
    subprocess.call([command], shell = True)

    # Move input and output files in the appropriate folder in SPheno
    if os.path.exists(SPHENO_work_folder + "Input_Files/" + SPheno_input_new_name):
        os.remove(SPHENO_work_folder + "Input_Files/" + SPheno_input_new_name)	
    shutil.move(SPHENO_folder + SPheno_input_new_name, SPHENO_work_folder + "Input_Files/")
    if not os.path.exists(SPHENO_work_folder + "Output_Files/"):
        os.mkdir(SPHENO_work_folder + "Output_Files/")
    if os.path.exists(SPHENO_work_folder + "Output_Files/" + SPheno_output_new_name):
        os.remove(SPHENO_work_folder + "Output_Files/" + SPheno_output_new_name)
    shutil.move(SPHENO_folder + SPheno_output_new_name, SPHENO_work_folder + "Output_Files/")

# Remove undesired files
output_list = ["BR_H_NP.dat", "BR_Hplus.dat", "BR_t.dat", "effC.dat", "LEP_HpHm_CS_ratios.dat", "Messages.out", "MH_GammaTot.dat", "MHplus_GammaTot.dat", "SPheno.out", model_name + "_low.out", "WC." + model_name + "_1.json", "WC." + model_name + "_2.json", "WHIZARD.par." + model_name]
for files in output_list:
    if os.path.exists(SPHENO_folder + files):
        os.remove(SPHENO_folder + files)

if os.path.exists(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name):
    os.remove(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name)

print("\n")
print("All SPheno spectra produced\n")
print("\n")


####################
### 3. RUN MICROMEGA
####################
## The micrOMEGA executable file is already available and compiled from previous script
print("Run micrOMEGA\n")
micrOMEGA_input_name = "SPheno.spc." + model_name

## Chose preferred micrOMEGA executable and output name. Use Simon to allign with his outputs
if Simon_micrOMEGA:
    micrOMEGA_exe_file = "micrOMEGA_" + model_name + "_run_Simon"
    micrOMEGA_output_name_default = "micromegas_data.out"
    micrOMEGA_output_channels_default = "micromegas_channels.out"
    micrOMEGA_output_spheno_default = "spheno_data.out"
else:
    micrOMEGA_exe_file = "micrOMEGA_" + model_name + "_run"
    micrOMEGA_output_name_default = "micrOMEGA_output.out"


command = "./"+micrOMEGA_exe_file
os.chdir(mircOMEGA_work_folder)
# Delete previous micrOMEGA output
if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_name_default):
	os.remove(mircOMEGA_work_folder + micrOMEGA_output_name_default)

for files in SPheno_output_list:
    # Remove existent input file
    if os.path.exists(mircOMEGA_work_folder + micrOMEGA_input_name):
        os.remove(mircOMEGA_work_folder + micrOMEGA_input_name)
    # Copy SPheno output as micrOMEGA input
    shutil.copy(SPHENO_work_folder + "Output_Files/" + files, mircOMEGA_work_folder + micrOMEGA_input_name)

    # Run micrOMEGA
    print("Running micrOMEGA for " + files + "\n")
    subprocess.call([command], shell = True)

print("micrOMEGA parameter scan complete\n")
micrOMEGA_output_name = "micrOMEGA_output_"
if len(parameter_names_scan) and len(block_names_scan):
    for i in range(0, len(parameter_names_scan)):
        micrOMEGA_output_name = micrOMEGA_output_name + parameter_names_scan[i] + "_" + str(parameter_table_scan[i][0]) + "_" + str(parameter_table_scan[0][-1]) + "_"
    for i in range(0, len(block_names_scan)):
        micrOMEGA_output_name = micrOMEGA_output_name + block_names_scan[i] + "_" + str(block_entries_scan[i][0][0][0]) + "_" + str(block_entries_scan[i][-1][0][0]) + "_"
elif len(parameter_names_scan):
    for i in range(0, len(parameter_names_scan)):
        micrOMEGA_output_name = micrOMEGA_output_name + parameter_names_scan[i] + "_" + str(parameter_table_scan[i][0]) + "_" + str(parameter_table_scan[0][-1]) + "_"
elif len(block_names_scan):
    for i in range(0, len(block_names_scan)):
        micrOMEGA_output_name = micrOMEGA_output_name + block_names_scan[i] + "_" + str(block_entries_scan[i][0][0][0]) + "_" + str(block_entries_scan[i][-1][0][0]) + "_"
micrOMEGA_output_name = micrOMEGA_output_name[:-1] + ".out"
if os.path.exists(local_path + micrOMEGA_output_name):
    os.remove(local_path + micrOMEGA_output_name)
shutil.move(mircOMEGA_work_folder + micrOMEGA_output_name_default, local_path + micrOMEGA_output_name)

if Simon_micrOMEGA:
    micrOMEGA_output_channels = "micrOMEGA_channels_"
    if len(parameter_names_scan) and len(block_names_scan):
        for i in range(0, len(parameter_names_scan)):
            micrOMEGA_output_channels = micrOMEGA_output_channels + parameter_names_scan[i] + "_" + str(parameter_table_scan[i][0]) + "_" + str(parameter_table_scan[0][-1]) + "_"
        for i in range(0, len(block_names_scan)):
            micrOMEGA_output_channels = micrOMEGA_output_channels + block_names_scan[i] + "_" + str(block_entries_scan[i][0][0][0]) + "_" + str(block_entries_scan[i][-1][0][0]) + "_"
    elif len(parameter_names_scan):
        for i in range(0, len(parameter_names_scan)):
            micrOMEGA_output_channels = micrOMEGA_output_channels + parameter_names_scan[i] + "_" + str(parameter_table_scan[i][0]) + "_" + str(parameter_table_scan[0][-1]) + "_"
    elif len(block_names_scan):
        for i in range(0, len(block_names_scan)):
            micrOMEGA_output_channels = micrOMEGA_output_channels + block_names_scan[i] + "_" + str(block_entries_scan[i][0][0][0]) + "_" + str(block_entries_scan[i][-1][0][0]) + "_"
    micrOMEGA_output_channels = micrOMEGA_output_channels[:-1] + ".out"
    if os.path.exists(local_path + micrOMEGA_output_channels):
        os.remove(local_path + micrOMEGA_output_channels)
    shutil.move(mircOMEGA_work_folder + micrOMEGA_output_channels_default, local_path + micrOMEGA_output_channels)

    micrOMEGA_output_spheno = "micrOMEGA_spheno_"
    if len(parameter_names_scan) and len(block_names_scan):
        for i in range(0, len(parameter_names_scan)):
            micrOMEGA_output_spheno = micrOMEGA_output_spheno + parameter_names_scan[i] + "_" + str(parameter_table_scan[i][0]) + "_" + str(parameter_table_scan[0][-1]) + "_"
        for i in range(0, len(block_names_scan)):
            micrOMEGA_output_spheno = micrOMEGA_output_spheno + block_names_scan[i] + "_" + str(block_entries_scan[i][0][0][0]) + "_" + str(block_entries_scan[i][-1][0][0]) + "_"
    elif len(parameter_names_scan):
        for i in range(0, len(parameter_names_scan)):
            micrOMEGA_output_spheno = micrOMEGA_output_spheno + parameter_names_scan[i] + "_" + str(parameter_table_scan[i][0]) + "_" + str(parameter_table_scan[0][-1]) + "_"
    elif len(block_names_scan):
        for i in range(0, len(block_names_scan)):
            micrOMEGA_output_spheno = micrOMEGA_output_spheno + block_names_scan[i] + "_" + str(block_entries_scan[i][0][0][0]) + "_" + str(block_entries_scan[i][-1][0][0]) + "_"
    micrOMEGA_output_spheno = micrOMEGA_output_spheno[:-1] + ".out"
    if os.path.exists(local_path + micrOMEGA_output_spheno):
        os.remove(local_path + micrOMEGA_output_spheno)
    shutil.move(mircOMEGA_work_folder + micrOMEGA_output_spheno_default, local_path + micrOMEGA_output_spheno)




#if Simon_micrOMEGA:
    #if os.path.exists(local_path + micrOMEGA_output_name_channels):
        #os.remove(local_path + micrOMEGA_output_name_channels)
    #shutil.move(mircOMEGA_work_folder + micrOMEGA_output_name_channels, local_path + micrOMEGA_output_name_channels)
    #if os.path.exists(local_path + micrOMEGA_output_name_spheno):
        #os.remove(local_path + micrOMEGA_output_name_spheno)
    #shutil.move(mircOMEGA_work_folder + micrOMEGA_output_name_spheno, local_path + micrOMEGA_output_name_spheno)

print("Results available in the " + micrOMEGA_output_name)
if Simon_micrOMEGA:
    print("Info on the channels available in " + micrOMEGA_output_channels)
    print("Info on the LFV observables available in " + micrOMEGA_output_spheno)

print("\n")

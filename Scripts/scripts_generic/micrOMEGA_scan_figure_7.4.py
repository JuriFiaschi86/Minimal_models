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
import numpy as np

### 0 -> Preliminary settings
# Setting model names
model_name = "T13Balpha0"
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


# Setting the fixed parameters
## Give the name and the values of the fixed input parameters
## The parameter names must match the SPheno input file names
SPheno_input_template_name = "LesHouches.in." + model_name + "_low"
SPheno_new_template_name = "Fixed_LesHouches.in." + model_name

parameter_names = ["lambdaInput","mphiInput","mCPsiInput","mpsipsipInput","lambda1Input","lambda3Input","lambda4Input","lambda5Input"]
parameter_table = [0.2612, 0, 1000000, 1000000, 0, 0, 0, 0]

block_name = ["LAM6IN"]
block_entries = [0, 0, 0]

SPheno_input_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_input_template_name, "r")
text = SPheno_input_template.read()
SPheno_input_template.close()
SPheno_input_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_input_template_name, "r")
SPheno_new_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name, "w")


old_parameter_line = []
new_parameter_line = []
parameter_value_string_old = []
parameter_value_string_new = []

command = "SPheno_new_template.write(line"
for i in range(0, len(parameter_names)):
    old_parameter_line.append(re.findall("^.*" + parameter_names[i] + ".*$", text, re.MULTILINE)[0])
    if int(old_parameter_line[i][:4]) < 10:
        space = 5
    else:
        space = 6
    parameter_value_string_old.append(old_parameter_line[i][space:][:13])
    parameter_value_string_new.append("{:.7E}".format(parameter_table[i]))
    new_parameter_line.append(old_parameter_line[i].replace(parameter_value_string_old[i], parameter_value_string_new[i]))
		
    command = command + ".replace(old_parameter_line["+ str(i) +"], new_parameter_line["+ str(i) +"])"
    
    
old_parameter_line_block = []
new_parameter_line_block = []
parameter_value_string_old_block = []
parameter_value_string_new_block = []

for j in range(0, len(block_entries)):
    parameter_name_block = block_name[0].lower()[:-2] + "\(" + str(j+1) + "\)"
    old_parameter_line_block.append(re.findall("^.*" + parameter_name_block + ".*$", text, re.MULTILINE)[0])
    parameter_value_string_old_block.append(old_parameter_line_block[j][5:][:13])
    parameter_value_string_new_block.append("{:.7E}".format(block_entries[j]))
    new_parameter_line_block.append(old_parameter_line_block[j].replace(parameter_value_string_old_block[j], parameter_value_string_new_block[j]))
    command = command + ".replace(old_parameter_line_block["+ str(j) +"], new_parameter_line_block["+ str(j) +"])"
    
command = command + ")"

for line in SPheno_input_template:
    exec(command)

SPheno_input_template.close()
SPheno_new_template.close()



# Settings of the parameter scan
## Give the name and the values of the parameters to scan over (as many as wanted)
## The parameter names must match the SPheno input file names

parameter_names = ["mphiInput","lambda1Input"]
parameter_table = [(100, 500, 1000, 1500, 2000, 3000, 4000, 5000),(0.1, 1)]

combine_param = list(itertools.product(*parameter_table))

#SPheno_input_template_name = "LesHouches.in." + model_name + "_low"
for i in range(0, len(combine_param)):
	# Preparing the input files
	SPheno_input_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name, "r")
	text = SPheno_input_template.read()
	SPheno_input_template.close()
	SPheno_input_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name, "r")

	SPheno_input_new_name = ""
	for j in range(0, len(parameter_names)):
		SPheno_input_new_name = SPheno_input_new_name + parameter_names[j] + "_" + str(combine_param[i][j]) + "_" 
	SPheno_input_new_name = SPheno_input_new_name[:-1] + ".in"
	SPheno_input_new = open(SPHENO_folder + SPheno_input_new_name, "w")

	# Changing the flags
	old_line_flag_50 = " 50 1               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) "
	old_line_flag_77 = " 77 0               # Output for MicrOmegas (running masses for light quarks; real mixing matrices)   "
	new_line_flag_50 = " 50 0               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) "
	new_line_flag_77 = " 77 1               # Output for MicrOmegas (running masses for light quarks; real mixing matrices)   "
	
	# Changing the value of parameter of the scan
	old_parameter_line = []
	new_parameter_line = []
	parameter_value_string_old = []
	parameter_value_string_new = []
	command = "SPheno_input_new.write(line.replace(old_line_flag_50, new_line_flag_50).replace(old_line_flag_77, new_line_flag_77)"
	for j in range(0, len(parameter_names)):
		old_parameter_line.append(re.findall("^.*" + parameter_names[j] + ".*$", text, re.MULTILINE)[0])
		if int(old_parameter_line[j][:4]) < 10:
			space = 5
		else:
			space = 6
		parameter_value_string_old.append(old_parameter_line[j][space:][:13])
		parameter_value_string_new.append("{:.7E}".format(combine_param[i][j]))
		new_parameter_line.append(old_parameter_line[j].replace(parameter_value_string_old[j], parameter_value_string_new[j]))
		
		command = command + ".replace(old_parameter_line["+ str(j) +"], new_parameter_line["+ str(j) +"])"
	command = command + ")"
		
	for line in SPheno_input_template:
		exec(command)
	
	SPheno_input_template.close()
	SPheno_input_new.close()

	# Run SPheno for each input files
	print("\n")
	print("Producing SPheno spectrum for " + SPheno_input_new_name + "\n")
	SPheno_output_new_name = SPheno_input_new_name[:-3] + ".out"
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
output_list = ["BR_H_NP.dat","BR_Hplus.dat","BR_t.dat","effC.dat","LEP_HpHm_CS_ratios.dat","Messages.out","MH_GammaTot.dat","MHplus_GammaTot.dat","SPheno.out","T13Balpha0_low.out","WC.T13Balpha0_1.json","WC.T13Balpha0_2.json","WHIZARD.par.T13Balpha0"]
for files in output_list:
	if os.path.exists(SPHENO_folder + files):
		os.remove(SPHENO_folder + files)

if os.path.exists(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name):
    os.remove(SPHENO_work_folder + "Input_Files/" + SPheno_new_template_name)

print("\n")
print("All SPheno spectra produced\n")
print("\n")



### 1 -> Run micrOMEGA
## The micrOMEGA executable file is already available and compiled from previous script
print("Run micrOMEGA\n")
micrOMEGA_input_name = "SPheno.spc." + model_name
micrOMEGA_exe_file = "micrOMEGA_" + model_name + "_run"
command = "./"+micrOMEGA_exe_file
os.chdir(mircOMEGA_work_folder)
# Delete previous micrOMEGA output
micrOMEGA_output_name_default = "micrOMEGA_output.out"
if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_name_default):
	os.remove(mircOMEGA_work_folder + micrOMEGA_output_name_default)

for i in range(0, len(combine_param)):
	# Remove existent input file
	if os.path.exists(mircOMEGA_work_folder + micrOMEGA_input_name):
		os.remove(mircOMEGA_work_folder + micrOMEGA_input_name)
	# Copy SPheno output as micrOMEGA input
	SPheno_output_name = ""
	for j in range(0, len(parameter_names)):
		SPheno_output_name = SPheno_output_name + parameter_names[j] + "_" + str(combine_param[i][j]) + "_"
	SPheno_output_name = SPheno_output_name[:-1] + ".out"
	shutil.copy(SPHENO_work_folder + "Output_Files/" + SPheno_output_name, mircOMEGA_work_folder + micrOMEGA_input_name)

	# Run micrOMEGA
	print("Running micrOMEGA for " + SPheno_output_name + "\n")
	subprocess.call([command], shell = True)

print("micrOMEGA parameter scan complete\n")
micrOMEGA_output_name = "micrOMEGA_output_"
for j in range(0, len(parameter_names)):
	micrOMEGA_output_name = micrOMEGA_output_name + parameter_names[j] + "_" + str(parameter_table[j][0]) + "_" + str(parameter_table[j][-1]) + "_"
micrOMEGA_output_name = micrOMEGA_output_name[:-1] + ".out"
if os.path.exists(local_path + micrOMEGA_output_name):
	os.remove(local_path + micrOMEGA_output_name)
shutil.move(mircOMEGA_work_folder + micrOMEGA_output_name_default, local_path + micrOMEGA_output_name)
print("Results are available in the " + micrOMEGA_output_name + "\n")

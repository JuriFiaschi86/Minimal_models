#! /usr/bin/python3

# This script sets up the necessary files for a parameter scan with micrOMEGA.
# Run SARAH to generate SPheno source code and model files
# The SPheno source code generated my SARAH is copied in the SPheno folder and compiled
### An example of SPheno spectrum is generated
# A new project for micrOMEGA is initialized and the model files from SARAH are copied in the appropriate folder
# The personalised micrOMEGA code is compiled and ready to be called by another script for the parameter scan


import shutil
import subprocess
import os
import stat
import glob

# Setting model names
model_name_long = "T1-3-B_0_scalar2g"
model_name = "T13Balpha02g"

# Setting codes versions
SARAH_version = "4.13.0"
SPheno_version = "4.0.3"
micrOMEGA_version = "4.3.5"
Mathematica_version = "11.0.1"

# Setting main directories
SARAH_folder = "/local0/fiaschi/SARAH-" + SARAH_version
SPHENO_folder = "/local0/fiaschi/SPheno/"
mircOMEGA_folder = "/local0/fiaschi/micrOMEGAS/"
MathKernel_script = "/Applic.ZIV/Mathematica/" + Mathematica_version + "/Executables/MathKernel"
local_path = os.path.dirname(os.path.abspath(__file__)) + "/"

### 1 -> Run SARAH to generate model files and SPheno code
print()
print('================ RUNNING SARAH =================')
print()
print("Run SARAH to generate SPheno source code and model files\n")
os.chdir(local_path)
SARAH_file_name = "run_SARAH.m"
command = MathKernel_script + " -noprompt -run \"<<"+SARAH_file_name+"\""
subprocess.call([command], shell = True)

### 2 -> copy the output of SARAH into SPHENO folder
# check if the directory already exist, and in case delete it
SARAH_SPheno_output = SARAH_folder + "/Output/" + model_name_long + "/EWSB/SPheno/"
SPheno_model_path = SPHENO_folder + "SPheno-" + SPheno_version + "/" + model_name
if os.path.exists(SPheno_model_path):
	shutil.rmtree(SPheno_model_path)
# copy the directory and change name (using subprocess method)
print("Copying SARAH output folder into SPheno\n")
command = "cp -r " + SARAH_SPheno_output + " " + SPheno_model_path
subprocess.call([command], shell = True)
print("Folder copied into SPheno\n")


### 3 -> compile the SPheno code
print()
print('=============== COMPILING SPHENO ===============')
print()
print("Compiling the SPheno code\n")
os.chdir(SPHENO_folder + "SPheno-"+ SPheno_version)
command = "make Model="+model_name
subprocess.call([command], shell = True)
print("SPheno code compiled\n")


### 4 -> set up and compile micrOMEGA
print()
print('============= COMPILING MICROMEGAS =============')
print()
print("Creating micrOMEGA inputs\n")
micrOMEGA_path = mircOMEGA_folder + "micromegas_" + micrOMEGA_version +"/"
# Check if micrOMEGA new project folder already exist
if os.path.exists(micrOMEGA_path + model_name):
	print("micrOMEGA project already existent\n")
# Create micrOMEGA new project
else:
	os.chdir(micrOMEGA_path)
	command = "./newProject " + model_name
	subprocess.call([command], shell = True)
# Copy relevant files to micrOMEGA new project folder
print("Copying relevant files to micrOMEGA new project folder\n")
# Copy from SARAH
model_files = ["func1.mdl","lgrng1.mdl","prtcls1.mdl","vars1.mdl"]
CalcOmega_file = "micrOMEGA_" + model_name + "_run_Simon.c"
SARAH_model_output = SARAH_folder + "/Output/" + model_name_long + "/EWSB/CHep/"
micrOMEGA_model_path = micrOMEGA_path + model_name + "/work/models/"
# Check if files already exist, then delete it and put the new ones
for files in model_files:
	if os.path.exists(micrOMEGA_model_path + files):
    		os.remove(micrOMEGA_model_path + files)
	shutil.copy(SARAH_model_output + files, micrOMEGA_model_path)
if os.path.exists(micrOMEGA_model_path + CalcOmega_file):
    	os.remove(micrOMEGA_model_path + CalcOmega_file)
shutil.copy(local_path + CalcOmega_file, micrOMEGA_path + model_name)
# Compile mircOMEGA
print("Compile micrOMEGA\n")
os.chdir(micrOMEGA_path + model_name)
command = "make main="+CalcOmega_file
subprocess.call([command], shell = True)
print("micrOMEGA compiled successfully\n")
print("Ready to run the script for parameter scan\n")

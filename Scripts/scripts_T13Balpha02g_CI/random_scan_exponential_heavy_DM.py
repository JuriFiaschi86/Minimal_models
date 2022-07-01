#! /usr/bin/python3

# WRITE COMMENTS

import shutil
import subprocess
import os
import re
import math
import random
import numpy as np 

import Casas_Ibarra_function

########################
### 0. SETTING THE FLAGS
########################
## Select how many points do we want
max_success = 10000
max_iteractions = 1000000000

## Flags to define the purpose of the code

# Flag for Casas-Ibarra (0 = off / 1 = on)
Casas_Ibarra_flag = 1

# Flag for putting the same value in the diagonals of "MPHI2IN" (0 = different masses / 1 = same masses)
flag_equal_scalars_masses = 0

# Flag for putting equal values to all block Lambda1 components (0 = different couplings / 1 = same couplings)
flag_equal_couplings = 1

# Flag for scalar DM (0 = off / 1 = on)
flag_scalar_DM = 1

# Flag for LFV check (0 = off / 1 = on)
flag_LFV_check = 0
# Experimental constrains
mu_e_gamma_bound = 5.7E-13
tau_e_gamma_bound = 3.3E-08
tau_mu_gamma_bound = 4.4E-08
mu_3e_bound = 1.0E-12
tau_3mu_bound = 2.1E-08
tau_e_2mu_bound = 2.7E-08
tau_mu_2e_bound = 1.8E-08
tau_3e_bound = 2.7E-08
mu_Ti_bound = 4.3E-12
mu_Au_bound = 7.0E-13

# Flag for DM relic density check (0 = off / 1 = on)
flag_relic_density_check = 0
# Experimental constrains
measured_relic = 0.12
tollerance_relic = 0.01

# Flag for heavy scalar DM
flag_heavy_DM = 1
threshold_DM_mass = 1000


###########################
### 1. PRELIMINARY SETTINGS
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
results_folder = local_path + "Heavy_DM_results/"


####################################
### 2. DEFINE MODEL INPUT PARAMETERS
####################################

## Identify the model parameters
#model_parameters = ["lambdaInput", "mCPsiInput", "mpsipsipInput", "lambda4Input", "lambda5Input"]
#model_blocks = ["MPHI2IN", "LAM1IN", "LAM6IN"]

## Setting the fixed parameters
parameter_names_fixed = ["mCPsiInput","mpsipsipInput","lambda4Input","lambda5Input"]
parameter_table_fixed = [1.0E+03, 1.5E+03, 5.0E-03, 5.0E-03]

block_names_fixed = ["MPHI2IN", "LAM1IN", "LAM6IN"]
block_entries_fixed = [(((1.4E+03)**2, 0), (0, (3.0E+03)**2)), ((0, 0), (0, 0)), ((0.0003, 0.0045), (0.0003, 0.0045), (0.003, 0.0045))]

## Setting the parameters of the scan

# Setting the name of the parameters of the scan
parameter_names_scan = ["mCPsiInput","mpsipsipInput","lambda4Input","lambda5Input"]

# Setting the range of the parameters of the scan (min/max of each parameter)
parameter_table_scan = [(1, 4), (1, 4), (-6, 0), (-6, 0)]

#parameter_names_scan = []
#parameter_table_scan = []


## Setting the blocks of the scan

# Setting the name of the blocks of the scan
block_names_scan = ["MPHI2IN", "LAM1IN"]

# Setting the range of the blocks of the scan (min/max of each block)
# (for mphi2 I set the limits for the value of the mass, then it gets squared afterwards)
# (also the limits are set such that ((mphi_11/22_min, mphi_11/22_max), (mphi_12/21_min, mphi_12/21_max)))
block_entries_scan = [((1, 4), (-6, 4)), (-6, 0)]

#block_names_scan = []
#block_entries_scan = []


# Number of generations of scalars (> 1 necessary for neutrino masses)
n_generation = 2

# Name of the necessary parameters for Casas-Ibarra parametrisation
Casas_Ibarra_input = ["mCPsiInput", "mpsipsipInput", "MPHI2IN", "LAM1IN", "lambda4Input", "lambda5Input"]
Casas_Ibarra_output = "LAM6IN"
Casas_Ibarra_input_values = []


##############################################
### 3. DEFINE FUNCTIONS TO WRITE SPHENO INPUTS
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

## Function to write random parameter value (masses linear random, couplings exponential random)
#def write_spheno_random_parameter (input_text, parameter_name, parameter_table):
    #if parameter_name[0] == "m":
        #random_value = random.uniform(parameter_table[0], parameter_table[1])
    #elif parameter_name[0:3] == "lam":
        #random_value = math.pow(10, random.uniform(parameter_table[0], parameter_table[1]))
    #input_text = write_spheno_input_param(input_text, parameter_name, random_value)
    #return input_text
    
## Function to write random parameter value (both masses and couplings exponential random)
#def write_spheno_random_parameter (input_text, parameter_name, parameter_table):
    #random_value = math.pow(10, random.uniform(parameter_table[0], parameter_table[1]))
    #input_text = write_spheno_input_param(input_text, parameter_name, random_value)
    #return input_text
    
# Function to write random parameter value (both masses and couplings exponential random, lambda 1 and 5 can be negative)
def write_spheno_random_parameter (input_text, parameter_name, parameter_table):
    if parameter_name == "lambda5Input" or parameter_name == "lambdaInput":
        random_value = ((-1)**random.randint(0,1))*math.pow(10, random.uniform(parameter_table[0], parameter_table[1]))
    else:
        random_value = math.pow(10, random.uniform(parameter_table[0], parameter_table[1]))
    input_text = write_spheno_input_param(input_text, parameter_name, random_value)
    return input_text


# Function to write random block entries
def write_spheno_random_block (input_text, block_name, block_entry):
    if block_name == "MPHI2IN":
        if flag_equal_scalars_masses:
            #random_value = random.uniform(block_entry[0], block_entry[1])
            random_value = math.pow(math.pow(10, random.uniform(block_entry[0][0], block_entry[0][1])), 2)
            random_block = (((random_value, 0), (0, random_value)))
        else:
            random_value_1 = math.pow(math.pow(10, random.uniform(block_entry[0][0], block_entry[0][1])), 2)
            random_value_2 = math.pow(math.pow(10, random.uniform(block_entry[0][0], block_entry[0][1])), 2)
            random_value_off = math.pow(math.pow(10, random.uniform(block_entry[1][0], block_entry[1][1])), 2)
            random_block = ((random_value_1, random_value_off), (random_value_off, random_value_2))
    elif block_name == "LAM1IN":
        if flag_equal_couplings:
            random_value = math.pow(10, random.uniform(block_entry[0], block_entry[1]))
            random_block = ((random_value, random_value), (random_value, random_value))
        else:
            random_value = math.pow(10, random.uniform(block_entry[0], block_entry[1]))
            random_block = ((random_value * random.uniform(0, 1), random_value * random.uniform(0, 1)), (random_value * random.uniform(0, 1), random_value * random.uniform(0, 1)))
    
    input_text = write_spheno_input_block(input_text, block_name, random_block)
    return input_text

# Function to change the flags in the Spheno input file for micrOMEGA run
def write_spheno_flag(input_text):
    old_line_flag_13 = " 13 1               # 3-Body decays: none (0), fermion (1), scalar (2), both (3) "
    new_line_flag_13 = " 13 0               # 3-Body decays: none (0), fermion (1), scalar (2), both (3) "
    old_line_flag_16 = " 16 1               # One-loop decays "
    new_line_flag_16 = " 16 0               # One-loop decays "
    old_line_flag_50 = " 50 1               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) "
    new_line_flag_50 = " 50 0               # Majorana phases: use only positive masses (put 0 to use file with CalcHep/Micromegas!) "
    old_line_flag_55 = " 55 0               # Calculate loop corrected masses "
    new_line_flag_55 = " 55 1               # Calculate loop corrected masses "
    #old_line_flag_67 = " 67 1               # effective Higgs mass calculation "
    #new_line_flag_67 = " 67 0               # effective Higgs mass calculation "
    old_line_flag_77 = " 77 0               # Output for MicrOmegas (running masses for light quarks; real mixing matrices)   "
    new_line_flag_77 = " 77 1               # Output for MicrOmegas (running masses for light quarks; real mixing matrices)   "
    input_text = input_text.replace(old_line_flag_13, new_line_flag_13)
    input_text = input_text.replace(old_line_flag_16, new_line_flag_16)
    input_text = input_text.replace(old_line_flag_50, new_line_flag_50)
    input_text = input_text.replace(old_line_flag_55, new_line_flag_55)
    #input_text = input_text.replace(old_line_flag_67, new_line_flag_67)
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


# Function to check if the point has scalar DM. Reads SPheno output and returns True if the lightest neutral particle is scalar, False if the lightest neutral particle is fermion
def check_scalar_DM(input_text):
    m_eta0_1 = float([line for line in input_text.split('\n') if "# eta0_1" in line][0][15:][:14])
    m_eta0_2 = float([line for line in input_text.split('\n') if "# eta0_2" in line][0][15:][:14])
    m_chi_1 = float([line for line in input_text.split('\n') if "# Fchi_1" in line][0][15:][:14])
    m_chi_2 = float([line for line in input_text.split('\n') if "# Fchi_2" in line][0][15:][:14])
    m_chi_3 = float([line for line in input_text.split('\n') if "# Fchi_3" in line][0][15:][:14])
    check = min(m_eta0_1, m_eta0_2) < min(m_chi_1, m_chi_2, m_chi_3)
    return check

# Function to select heavy DM. Reads SPheno output and returns True if lightest scalar mass is > 2 TeV, False if lightest scalar mass is < 2 TeV
def check_heavy_DM(input_text):
    m_eta0_1 = float([line for line in input_text.split('\n') if "# eta0_1" in line][0][15:][:14])
    m_eta0_2 = float([line for line in input_text.split('\n') if "# eta0_2" in line][0][15:][:14])
    check = min(m_eta0_1, m_eta0_2) > threshold_DM_mass
    return check


## Function to write headers of "micromegas_data" and "spheno_data"
#def write_headers():
    #header_micromegas_data = "# mϕ₁₁\tmϕ₁₂\tmϕ₂₁\tmϕ₂₂\tmΨ\tmψψ'\tλ₁₁₁\tλ₁₁₂\tλ₁₂₁\tλ₁₂₂\tλ₄\tλ₅\tλ₆₁₁\tλ₆₂₁\tλ₆₃₁\tλ₆₁₂\tλ₆₂₂\tλ₆₃₂\tλ\tv\tmd1\tmd2\tmd3\tmu1\tmu2\tmu3\tme1\tme2\tme3\tmν1\tmν2\tmν3\tmZ\tmW\tmH\tmη₁+\tmη₂+\tmη₁0\tmη₂0\tmψ\tmχ1\tmχ2\tmχ3\tΓZ\tΓW\tΓH\tΩh²\tXf\tσp(SI)\tσp(SD)"
    #header_spheno_data = "# mϕ₁₁\tmϕ₁₂\tmϕ₂₁\tmϕ₂₂\tmΨ\tmψψ'\tλ₁₁₁\tλ₁₁₂\tλ₁₂₁\tλ₁₂₂\tλ₄\tλ₅\tλ₆₁₁\tλ₆₂₁\tλ₆₃₁\tλ₆₁₂\tλ₆₂₂\tλ₆₃₂\tλ\tv\tBR(μ→eγ)\tBR(τ→eγ)\tBR(τ→μγ)\tBR(μ→3e)\tBR(τ→3e)\tBR(τ→3μ)"
    #with open(results_folder + micrOMEGA_output_data_final, "a") as micrOMEGA_data_final:
            #micrOMEGA_data_final.write(header_micromegas_data)
    #with open(results_folder + micrOMEGA_output_spheno_final, "a") as micrOMEGA_spheno_final:
            #micrOMEGA_spheno_final.write(header_spheno_data)


##################
### INITIALIZATION
##################

SPheno_input_temp_name = "temp.in"
SPheno_output_temp_name = SPheno_input_temp_name[:-3] + ".out"

error_line_1 = "Stopping calculation because of negative mass squared."
error_line_2 = "Mass spectrum converged, but negative mass squared present."
error_line_3 = "Error appeared in calculation of masses"    
error_line_4 = "in the calculation of the masses occurred a negative mass squared"
error_line_5 = "NaN appearing in MHp2"

# Define file names (uses Simon micrOMEGA source)
micrOMEGA_exe_file = "micrOMEGA_" + model_name + "_run_Simon"

micrOMEGA_output_data_default = "micromegas_data.out"
micrOMEGA_output_channels_default = "micromegas_channels.out"
micrOMEGA_output_spheno_default = "spheno_data.out"

micrOMEGA_output_data_final = "micromegas_data_final.out"
micrOMEGA_output_channels_final = "micromegas_channels_final.out"
micrOMEGA_output_spheno_final = "spheno_data_final.out"

### Wirte the headers of the results files
## Prepare a folder to store results
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
header_micromegas_data = "# mϕ₁₁\tmϕ₁₂\tmϕ₂₁\tmϕ₂₂\tmΨ\tmψψ'\tλ₁₁₁\tλ₁₁₂\tλ₁₂₁\tλ₁₂₂\tλ₄\tλ₅\tλ₆₁₁\tλ₆₂₁\tλ₆₃₁\tλ₆₁₂\tλ₆₂₂\tλ₆₃₂\tλ\tv\tmd1\tmd2\tmd3\tmu1\tmu2\tmu3\tme1\tme2\tme3\tmν1\tmν2\tmν3\tmZ\tmW\tmH\tmη₁+\tmη₂+\tmη₁0\tmη₂0\tmψ\tmχ1\tmχ2\tmχ3\tΓZ\tΓW\tΓH\tΩh²\tXf\tσp(SI)\tσp(SD)\n"
header_spheno_data = "# mϕ₁₁\tmϕ₁₂\tmϕ₂₁\tmϕ₂₂\tmΨ\tmψψ'\tλ₁₁₁\tλ₁₁₂\tλ₁₂₁\tλ₁₂₂\tλ₄\tλ₅\tλ₆₁₁\tλ₆₂₁\tλ₆₃₁\tλ₆₁₂\tλ₆₂₂\tλ₆₃₂\tλ\tv\tBR(μ→eγ)\tBR(τ→eγ)\tBR(τ→μγ)\tBR(μ→3e)\tBR(τ→3e)\tBR(τ→3μ)\n"
with open(results_folder + micrOMEGA_output_data_final, "a") as micrOMEGA_data_final:
        micrOMEGA_data_final.write(header_micromegas_data)
with open(results_folder + micrOMEGA_output_spheno_final, "a") as micrOMEGA_spheno_final:
        micrOMEGA_spheno_final.write(header_spheno_data)


###########################
### START OF THE ITERACTION
###########################


discarded_perturbativity = 0
discarded_LFV = 0
discarded_relic_density = 0
discarded_DM_candidate = 0
discarded_light_DM_candidate = 0
discarded_SPheano_fail = 0
discarded_micrOMEGA_fail = 0

n_success = 0
iteration = 1

### Loop conditions. I can enforce a limit on the sole number of positive results, or also on the number of interactios.

#while n_success < max_success and iteration <= max_iteractions:

while n_success < max_success:

###############################
### 4. WRITE SPHENO INPUT FILES
###############################

        
    print("### SUMMARY ###")
    print("Number of iterations = " + str(iteration-1))
    print("Accepted number of points = " + str(n_success))
    if flag_scalar_DM:
        print("Discarded for absence of scalar DM = " + str(discarded_DM_candidate))
    if flag_heavy_DM:
        print("Discarded for light scalar DM = " + str(discarded_light_DM_candidate))
    print("Discarded for perturbativity = " + str(discarded_perturbativity))
    if flag_LFV_check:
        print("Discarded for LFV = " + str(discarded_LFV))
    if flag_relic_density_check:
        print("Discarded for relic density = " + str(discarded_relic_density))
    print("Discarded for SPheno fail = " + str(discarded_SPheano_fail))
    print("Discarded for micrOMEGA fail = " + str(discarded_micrOMEGA_fail))
   
    print("\nGenerating new random point in parameter space.\n")

    # Open template Spheno input file
    SPheno_input_template = open(SPHENO_work_folder + "Input_Files/" + SPheno_input_template_name, "r")
    text = SPheno_input_template.read()
    SPheno_input_template.close()

    # Write the flags
    text = write_spheno_flag(text)
    # Write fixed parameters
    for i in range(0, len(parameter_names_fixed)):
        text = write_spheno_input_param(text, parameter_names_fixed[i], parameter_table_fixed[i])
    # Write fixed blocks
    for i in range(0, len(block_names_fixed)):
        text = write_spheno_input_block(text, block_names_fixed[i], block_entries_fixed[i])

    # Write random parameters
    for i in range(0, len(parameter_names_scan)):
        text = write_spheno_random_parameter(text, parameter_names_scan[i], parameter_table_scan[i])
    # Write random blocks
    for i in range(0, len(block_names_scan)):
        text = write_spheno_random_block(text, block_names_scan[i], block_entries_scan[i])
    
    
    ### Check for scalar DM (check that the scalar mass is smaller than fermion masses)
    if flag_scalar_DM:
        scalar_mass_1 = math.sqrt(float((re.findall("^.*mphi2\(1,1\).*$", text, re.MULTILINE)[0])[7:][:13]))
        scalar_mass_2 = math.sqrt(float((re.findall("^.*mphi2\(2,2\).*$", text, re.MULTILINE)[0])[7:][:13]))
        fermion_mass_1 = float((re.findall("^.*mCPsiInput.*$", text, re.MULTILINE)[0])[6:][:13])
        fermion_mass_2 = float((re.findall("^.*mpsipsipInput.*$", text, re.MULTILINE)[0])[6:][:13])
        
        if min(scalar_mass_1, scalar_mass_2) > min(fermion_mass_1, fermion_mass_2):
            print("No scalar DM candidate. Discard point.\n")
            discarded_DM_candidate += 1
            iteration += 1
            continue

    ## Write and check for Casas-Ibarra output
    if Casas_Ibarra_flag:
        Casas_Ibarra_output_value = write_Casas_Ibarra_parametrisation(text)
        text = write_spheno_input_block(text, Casas_Ibarra_output, Casas_Ibarra_output_value)
        if max(abs(Casas_Ibarra_output_value[0][0]),abs(Casas_Ibarra_output_value[0][1]),abs(Casas_Ibarra_output_value[1][0]),abs(Casas_Ibarra_output_value[1][1]),abs(Casas_Ibarra_output_value[2][0]),abs(Casas_Ibarra_output_value[2][1])) >= 4*math.pi:
            print("Lambda6 violates perturbativity. Discard point.\n")
            discarded_perturbativity += 1
            iteration += 1
            continue

    if os.path.exists(SPHENO_folder + SPheno_input_temp_name):
            os.remove(SPHENO_folder + SPheno_input_temp_name)

    SPheno_input_temp = open(SPHENO_folder + SPheno_input_temp_name, "w")
    SPheno_input_temp.write(text)
    SPheno_input_temp.close()
    

#################
### 5. RUN SPHENO
#################

    # Run SPheno
    if os.path.exists(SPHENO_folder + SPheno_output_temp_name):
            os.remove(SPHENO_folder + SPheno_output_temp_name)

    os.chdir(SPHENO_folder)
    command = SPHENO_folder + "SPheno-" + SPheno_version +"/bin/SPheno" + model_name + " " + SPHENO_folder + SPheno_input_temp_name + " " + SPHENO_folder + SPheno_output_temp_name
    
    SPheno_call = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (SPheno_output, SPheno_error) = SPheno_call.communicate()

    SPheno_call_status = SPheno_call.wait()

    SPheno_call_check_1 = (str(SPheno_output).find(error_line_1))
    SPheno_call_check_2 = (str(SPheno_output).find(error_line_2))
    SPheno_call_check_3 = (str(SPheno_output).find(error_line_3))
    SPheno_call_check_4 = (str(SPheno_output).find(error_line_4))
    SPheno_call_check_5 = (str(SPheno_output).find(error_line_5))

    if SPheno_call_check_1 != -1 or SPheno_call_check_2 != -1 or SPheno_call_check_3 != -1 or SPheno_call_check_4 != -1 or SPheno_call_check_5 != -1:
        print("SPheno output error. Discard point.\n")
        discarded_SPheano_fail += 1
        iteration += 1
        continue
    
    # Remove undesired files
    output_list = ["BR_H_NP.dat", "BR_Hplus.dat", "BR_t.dat", "effC.dat", "LEP_HpHm_CS_ratios.dat", "Messages.out", "MH_GammaTot.dat", "MHplus_GammaTot.dat", "SPheno.out", model_name + "_low.out", "WC." + model_name + "_1.json", "WC." + model_name + "_2.json", "WHIZARD.par." + model_name]
    for files in output_list:
        if os.path.exists(SPHENO_folder + files):
            os.remove(SPHENO_folder + files)
            
    
    ### Checks on DM mass
    if (flag_scalar_DM or flag_heavy_DM):
        SPheno_output = open(SPHENO_folder + SPheno_output_temp_name, "r")
        SPheno_output_text = SPheno_output.read()
        SPheno_output.close()
    
    ### Check for scalar DM
    if flag_scalar_DM and not(check_scalar_DM(SPheno_output_text)):
        print("No scalar DM candidate. Discard point.\n")
        discarded_DM_candidate += 1
        iteration += 1
        continue
    
    ### Check for heavy DM
    if flag_heavy_DM and not(check_heavy_DM(SPheno_output_text)):
        print("DM mass < " + str(threshold_DM_mass) + " GeV. Discard point.\n")
        discarded_light_DM_candidate += 1
        iteration += 1
        continue

        
###########################
### 6. CHECK LFV CONSTRAINS
###########################

    if flag_LFV_check:
        
        # Read LFV observable from SPheno    
        for line in open(SPHENO_folder + SPheno_output_temp_name, "r"):
            if "# BR(mu->e gamma)" in line:
                mu_e_gamma = float(line[12:26])
            if "# BR(tau->e gamma)" in line:
                tau_e_gamma = float(line[12:26])
            if "# BR(tau->mu gamma)" in line:
                tau_mu_gamma = float(line[12:26])
            if "# CR(mu-e, Ti)" in line:
                mu_Ti = float(line[12:26])
            if "# CR(mu-e, Au)" in line:
                mu_Au = float(line[12:26])
            if "# BR(mu->3e)" in line:
                mu_3e = float(line[12:26])
            if "# BR(tau->3e)" in line:
                tau_3e = float(line[12:26])
            if "# BR(tau->3mu)" in line:
                tau_3mu = float(line[12:26])
            if "# BR(tau- -> e- mu+ mu-)" in line:
                tau_e_2mu = float(line[12:26])
            if "# BR(tau- -> mu- e+ e-)" in line:
                tau_mu_2e = float(line[12:26])
                
        # Check of the constrains
        LFV_check = (mu_e_gamma <= mu_e_gamma_bound) and (tau_e_gamma <= tau_e_gamma_bound) and (tau_mu_gamma <= tau_mu_gamma_bound) and (mu_3e <= mu_3e_bound) and (tau_3mu <= tau_3mu_bound) and (tau_e_2mu <= tau_e_2mu_bound) and (tau_mu_2e <= tau_mu_2e_bound) and (tau_3e <= tau_3e_bound) and (mu_Ti <= mu_Ti_bound) and (mu_Au <= mu_Au_bound)
        if not LFV_check:
            print("LFV constrains failed. Discard point.\n")
            discarded_LFV += 1
            iteration += 1
            continue
            

####################
### 7. RUN MICROMEGA
####################

    micrOMEGA_input_name = "SPheno.spc." + model_name
    if os.path.exists(mircOMEGA_work_folder + micrOMEGA_input_name):
        os.remove(mircOMEGA_work_folder + micrOMEGA_input_name)
    # Copy SPheno output as micrOMEGA input
    shutil.copy(SPHENO_folder + SPheno_output_temp_name, mircOMEGA_work_folder + micrOMEGA_input_name)
    
    ### Prepare a folder to store results
    #if not os.path.exists(results_folder):
        #os.makedirs(results_folder)

    # Run micrOMEGA
    os.chdir(mircOMEGA_work_folder)
    command = "./"+micrOMEGA_exe_file
       
    micrOMEGA_call = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (micrOMEGA_output, micrOMEGA_error) = micrOMEGA_call.communicate()
    micrOMEGA_call_status = micrOMEGA_call.wait()
    if micrOMEGA_call_status:
        print("micrOMEGA output error. Discard point.\n")
        discarded_micrOMEGA_fail += 1
        iteration += 1
        continue
    
    ### Check relic density output
    if flag_relic_density_check:
        relic_density_upper_limit = measured_relic + tollerance_relic
        relic_density_lower_limit = measured_relic - tollerance_relic
        
        micrOMEGA_output = np.loadtxt(mircOMEGA_work_folder + micrOMEGA_output_data_default, delimiter="\t")
        
        relic_density = micrOMEGA_output[46]           
            
        if (relic_density < relic_density_upper_limit and relic_density > relic_density_lower_limit):
            print("Found a viable point at iteraction # " + str(iteration) + "\n")
            
            ### Store the positive result in the final files
            with open(mircOMEGA_work_folder + micrOMEGA_output_data_default, "r") as micrOMEGA_data_default:
                data_text = micrOMEGA_data_default.read()
            with open(results_folder + micrOMEGA_output_data_final, "a") as micrOMEGA_data_final:
                micrOMEGA_data_final.write(data_text)
            
            with open(mircOMEGA_work_folder + micrOMEGA_output_channels_default, "r") as micrOMEGA_channels_default:
                channels_text = micrOMEGA_channels_default.read()
            with open(results_folder + micrOMEGA_output_channels_final, "a") as micrOMEGA_channels_final:
                micrOMEGA_channels_final.write(channels_text)
                
            with open(mircOMEGA_work_folder + micrOMEGA_output_spheno_default, "r") as micrOMEGA_spheno_default:
                spheno_text = micrOMEGA_spheno_default.read()
            with open(results_folder + micrOMEGA_output_spheno_final, "a") as micrOMEGA_spheno_final:
                micrOMEGA_spheno_final.write(spheno_text)
            
            n_success += 1
            
        else:
            print("Relic density outside the bounds. Discard point.\n")            
            discarded_relic_density += 1
            
        ## Remove old files
        if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_data_default):
            os.remove(mircOMEGA_work_folder + micrOMEGA_output_data_default)
            
        if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_channels_default):
            os.remove(mircOMEGA_work_folder + micrOMEGA_output_channels_default)
            
        if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_spheno_default):
            os.remove(mircOMEGA_work_folder + micrOMEGA_output_spheno_default)
            
    else:
        ### Store the result in the final files
        with open(mircOMEGA_work_folder + micrOMEGA_output_data_default, "r") as micrOMEGA_data_default:
            data_text = micrOMEGA_data_default.read()
        with open(results_folder + micrOMEGA_output_data_final, "a") as micrOMEGA_data_final:
            micrOMEGA_data_final.write(data_text)
        
        with open(mircOMEGA_work_folder + micrOMEGA_output_channels_default, "r") as micrOMEGA_channels_default:
            channels_text = micrOMEGA_channels_default.read()
        with open(results_folder + micrOMEGA_output_channels_final, "a") as micrOMEGA_channels_final:
            micrOMEGA_channels_final.write(channels_text)
            
        with open(mircOMEGA_work_folder + micrOMEGA_output_spheno_default, "r") as micrOMEGA_spheno_default:
            spheno_text = micrOMEGA_spheno_default.read()
        with open(results_folder + micrOMEGA_output_spheno_final, "a") as micrOMEGA_spheno_final:
            micrOMEGA_spheno_final.write(spheno_text)
        
        n_success += 1
        
        ## Remove old files
        if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_data_default):
            os.remove(mircOMEGA_work_folder + micrOMEGA_output_data_default)
            
        if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_channels_default):
            os.remove(mircOMEGA_work_folder + micrOMEGA_output_channels_default)
            
        if os.path.exists(mircOMEGA_work_folder + micrOMEGA_output_spheno_default):
            os.remove(mircOMEGA_work_folder + micrOMEGA_output_spheno_default)
            
            
    if os.path.exists(mircOMEGA_work_folder + micrOMEGA_input_name):
        os.remove(mircOMEGA_work_folder + micrOMEGA_input_name)
    
    
###########################
### CONCLUDE THE ITERACTION
###########################

    iteration += 1
    

##################################
### 8. PRINT AND STORE THE RESULTS
##################################

print("\n")
print("### FINAL SUMMARY ###")
print("Number of iterations = " + str(iteration))
print("Accepted number of points = " + str(n_success))
if flag_scalar_DM:
    print("Discarded for absence of scalar DM = " + str(discarded_DM_candidate))
if flag_heavy_DM:
    print("Discarded for light scalar DM = " + str(discarded_light_DM_candidate))
print("Discarded for perturbativity = " + str(discarded_perturbativity))
if flag_LFV_check:
    print("Discarded for LFV = " + str(discarded_LFV))
if flag_relic_density_check:
    print("Discarded for relic density = " + str(discarded_relic_density))
print("Discarded for SPheno fail = " + str(discarded_SPheano_fail))
print("Discarded for micrOMEGA fail = " + str(discarded_micrOMEGA_fail))

print("\n")
print("### STORED RESULTS ###")
#print("SPheno outputs available in \"" + micrOMEGA_output_spheno_final + "\"")
#print("List of contributing channels for each data point available in \"" + micrOMEGA_output_channels_final + "\"")
#print("Data points available in \"" + micrOMEGA_output_data_final + "\"")
print("Files with stored results saved in \"" + results_folder + "\"")

#! /usr/bin/python3

# WRITE COMMENTS

import shutil
import subprocess
import os
import re
import random


print("Attempt to run SPheno")

path = "/home/users/fiaschi/Dokumente/SPheno"
os.chdir(path)

#command = "./SPheno-4.0.3/bin/SPhenoT13Balpha02g ./temp_error.in ./temp.out"
#command = "./SPheno-4.0.3/bin/SPhenoT13Balpha02g ./temp_fine.in ./temp.out"

command = "./SPheno-4.0.3/bin/SPhenoT13Balpha02g ./temp.in ./temp.out"


#subprocess.call([command], shell = True)


SPheno_call = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)


(SPheno_output, SPheno_error) = SPheno_call.communicate()

SPheno_status = SPheno_call.wait()
print("Command output : ", SPheno_output)
print("Command error : ", SPheno_error)
print("Command exit status/return code : ", SPheno_status)

error_line_1 = "Stopping calculation because of negative mass squared."
error_line_2 = "Mass spectrum converged, but negative mass squared present."
error_line_3 = "Error appeared in calculation of masses"
error_line_4 = "in the calculation of the masses occurred a negative mass squared"
error_line_5 = "NaN appearing in MHp2"
SPheno_call_check_1 = (str(SPheno_output).find(error_line_1))
SPheno_call_check_2 = (str(SPheno_output).find(error_line_2))
SPheno_call_check_3 = (str(SPheno_output).find(error_line_3))
SPheno_call_check_4 = (str(SPheno_output).find(error_line_4))
SPheno_call_check_5 = (str(SPheno_output).find(error_line_5))


print(SPheno_call_check_1)
print(SPheno_call_check_2)
print(SPheno_call_check_3)
print(SPheno_call_check_4)
print(SPheno_call_check_5)


print(SPheno_call_check_1 != -1 or SPheno_call_check_2 != -1 or SPheno_call_check_3 != -1 or SPheno_call_check_4 != -1 or SPheno_call_check_5 != -1)



#random_sign = random.choice('+-')
#print(random_sign)


#error_line = "Stopping calculation because of negative mass squared."
#print(str(output).find(error_line))

#error_line_2 = "Check"
#print(str(output).find(error_line_2))

##print(output.find(error_line))



#if not p_status:
    #print("Skipped")




#path = "/home/users/fiaschi/Dokumente/micrOMEGA/micromegas_4.3.5/T13Balpha02g"
#os.chdir(path)

#command = "./micrOMEGA_T13Balpha02g_run_Simon"

## Run micrOMEGA
##subprocess.call([command], shell = True)



#micrOMEGA_call = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
#(micrOMEGA_output, micrOMEGA_error) = micrOMEGA_call.communicate()

#micrOMEGA_call_status = micrOMEGA_call.wait()
#print("Command output : ", micrOMEGA_output)
#print("Command exit status/return code : ", micrOMEGA_call_status)

#if not micrOMEGA_call_status:
    #print("False")


#error_line = "Stopping calculation because of negative mass squared."
#SPheno_call_check = (str(SPheno_output).find(error_line))

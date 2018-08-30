#!/usr/bin/env python

# Write Defintions File for xx20

#----------------------------------------------------------------------
# Imports
#----------------------------------------------------------------------

import os
import sys
import argparse

#----------------------------------------------------------------------
# Class Defintions
#----------------------------------------------------------------------

class definition: 
    def __init__(self,scrit=0.005,load=0.01,charlen=0.0005,resolution=1e5):
        # Setting Default Variables
        self._scrit = scrit
        self._load  = load
        self._charlen = charlen
        self._resolution = resolution

    # Access Function
    def set_scrit(self,value):
	    self._scrit = value

    # Return functions
    def get_scrit(self):
        return str(self._scrit)
    def get_load(self):
        return str(self._load)
    def get_charlen(self):
        return str(self._charlen)
    def get_resolution(self):
        return str(self._resolution)

#----------------------------------------------------------------------
# Function Defintions
#----------------------------------------------------------------------
def myParser():
    """
    function : Checking Command line Arguments
    params   : null
    :return   : arguments
    """
    parser = argparse.ArgumentParser(description="Writes .def file from .dat file")

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-o', metavar='--Output', help=".def file out",dest='output',required=True)
    required.add_argument('-i', metavar='--Input', help=".dat file in",dest='input',required=True)

    parser.add_argument('-v', metavar='--Verbose', help="Increase Verbosity")
    args = parser.parse_args()

    return args



def read_dat(fname="xx20.dat"):
    dat_list = []
    dat_dict = {}

    with open(fname) as dat_file:
        # Get file infdrmation
        information = dat_file.readlines()

        # Clean information
        for i in information:
            line = i.split()
            for j in line:
                dat_list.append(j)
    
    # Add to dictionary for ease of use
    dat_dict.update({'element':dat_list[0]})
    dat_dict.update({'mesh':dat_list[1]})
    dat_dict.update({'partitioner':dat_list[2]})
    dat_dict.update({'nels':dat_list[3]})
    dat_dict.update({'nn':dat_list[4]})
    dat_dict.update({'nr':dat_list[5]})
    dat_dict.update({'nip':dat_list[6]})
    dat_dict.update({'nod':dat_list[7]})
    dat_dict.update({'loaded_nodes':dat_list[8]})
    dat_dict.update({'nres':dat_list[9]})
    dat_dict.update({'rho':dat_list[10]})
    dat_dict.update({'e':dat_list[11]})
    dat_dict.update({'v':dat_list[12]})
    dat_dict.update({'alpha1':dat_list[13]})
    dat_dict.update({'beta1':dat_list[14]})
    dat_dict.update({'nstep':dat_list[15]})
    dat_dict.update({'npri':dat_list[16]})
    dat_dict.update({'theta':dat_list[17]})
    dat_dict.update({'omega':dat_list[18]})
    dat_dict.update({'tol':dat_list[19]})
    dat_dict.update({'limit':dat_list[20]})

    return dat_dict

def writeToFile(defObj,datDict,fname):
    
    file = open(fname,"w")
 
    file.write("# Definitions File\n")
    file.write("\n")
    file.write(datDict["element"]+"\t#-Element\n")
    file.write(datDict["mesh"]+"\t#-Mesh\n")
    file.write(datDict["partitioner"]+"\t#-Partitioner\n")
    file.write(datDict["nels"]+"\t#-nels\n")
    file.write(datDict["nn"]+"\t#-nn\n")
    file.write(datDict["nr"]+"\t#-nr\n")
    file.write(datDict["nip"]+"\t#-nip\n")
    file.write(datDict["nod"]+"\t#-nod\n")
    file.write(datDict["loaded_nodes"]+"\t#-LoadedNodes\n")
    file.write(datDict["nres"]+"\t#-nres\n")
    file.write("\n")

    file.write("0.787400E+04\t#-Density\n")
    file.write("200.00000E+09\t#-YoungsModulus\n")
    file.write("0.400000E+00\t#-Poissons\n")
    file.write("0.000000E+00\t#-Alpha1\n")
    file.write("0.000000E+00\t#-Beta1\n")
    file.write("10\t\t#-nstep\n")
    file.write("1000\t\t#-npri\n")
    file.write("0.500000E+00\t#-theta\n")
    file.write("0.100000E-04\t#-Omega\n")
    file.write("0.100000E-05\t#-tol\n")
    file.write("10000\t\t#-limit\n")
    file.write("0.00001\t\t#-dtim\n")
    file.write("\n")
    file.write(defObj.get_charlen()+"\t#-characteristicLength\n")
    file.write(defObj.get_resolution()+"\t#-resolution\n")
    file.write(defObj.get_load()+"\t\t#-loadScale\n")
    file.write(defObj.get_scrit()+"\t\t#-scritScale\n")
    file.write("5\t\t#-crackStart\n")
    file.write("1.0\t\t#-exitDisplacement\n")

    file.close() 

#----------------------------------------------------------------------
# Main: Start Here
#----------------------------------------------------------------------

if __name__ == "__main__":
   
    print "Writing Defintions File"

    # Argument Parser
    args = myParser()

    if args.input is None:
        print("\nInput not set")
        print("Please see --help")
        exit()

    if args.output is None:
        print("\nOutput not set")
        print("Please see --help")
        exit()

    print("Input File Name  : {0}".format(args.input))
    print("Output File Name : {0}".format(args.output))

    # Read dat file to obtain mesh statistics
    datfile = read_dat(args.input)

    # Get Characteristic length
    # 1.25 * fe element size
    with open('xx20.mg','r') as f:
        data=f.read()
        data = data.split(' ') 
    
    char = float(data[5])*1.25
    dirname, filename = os.path.split(os.path.abspath(__file__)) 
    dirname=dirname.split('/')
    nodes = float(dirname[-1])
    
    if nodes == 0:
        cores = 1.0
    else:
        cores = nodes *24

    # Calculate Resolution
    # --------------------
    CELLS_CORE = 400000.0
    GRAINS = 512.0
    # --------------------
    total_cells = CELLS_CORE*cores
    res = int(total_cells/512.0)

    # Declare Definitions Object
    myDef = definition(scrit=1e-16,load=5000,charlen=char,resolution=res)                
    
    # write defintion to File
    writeToFile(myDef,datfile,args.output)







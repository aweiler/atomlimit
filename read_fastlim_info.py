import os
import subprocess
import sys
import logging
from math import *

#from collections import *
import yaml

from os import listdir
from os.path import  isdir,isfile, join

from collections import OrderedDict

class UnsortableList(list):
    def sort(self, *args, **kwargs):
        pass

class UnsortableOrderedDict(OrderedDict):
    def items(self, *args, **kwargs):
        return UnsortableList(OrderedDict.items(self, *args, **kwargs))

yaml.add_representer(UnsortableOrderedDict, yaml.representer.SafeRepresenter.represent_dict)




mypath = "/Users/aweiler/MC/AtomStatYaml/8TeV"

anaList = [f for f in listdir(mypath) if isdir(join(mypath, f))]

print anaList

for anaName in anaList:
    filename = "8TeV/" + anaName + "/SR_info.txt"
    fileA = open(filename, 'r')

    srlistSTR = []
    for lineIn in fileA:
        print lineIn
        srlistSTR.append( lineIn.split()  ) 
    del srlistSTR[0] # remove text header
    srlistSTR = [x for x in srlistSTR if x != []]  # remove empty
   
    print srlistSTR

    dictList = []

    for srlistEntry in srlistSTR:
        print srlistEntry

        sr_info = ''
        for idm in range(10, len(srlistEntry)): 
            sr_info += srlistEntry[idm]
            if not idm+1 == len(srlistEntry): sr_info +=' '
        print sr_info
        oneSRdict = dict()
        oneSRdict[sr_info] = dict (
            # Name = sr_info,
            Type = "CutAndCount",
            Luminosity =  float(srlistEntry[1]),
            rootS =  float(srlistEntry[0]),
            Observed =  int(srlistEntry[2]),
            Background =  float(srlistEntry[3]),
            SystematicError =  float(srlistEntry[4]),
            UpperLimit95obs =  float(srlistEntry[5]),
            UpperLimit95exp =  float(srlistEntry[6]),
            UpperXsecFB95obs =  float(srlistEntry[7]),
            UpperXsecFB95exp =  float(srlistEntry[8])
            ) 
        
        dictList.append(oneSRdict)

    yamlFile = dict( 
        Name= anaName, 
        SignalRegions = dictList
        )
    yamlName = anaName + ".stat"
    with open(yamlName, 'w') as outfile:
        outfile.write( yaml.dump(yamlFile,default_flow_style=False) )

    fileA.close()
    
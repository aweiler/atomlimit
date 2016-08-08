#! /usr/bin/env python

"""
"""

__author__ = "A.Weiler"
__version__ = "0.1"

# Add Control Region information ? 
#-> Check if flag


import os
import subprocess
import sys
import logging
from math import *

#from collections import *
import yaml
import json

from tabulate import tabulate

from os import listdir
from os.path import  isdir,isfile, join

from collections import OrderedDict

import operator

# Load
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

import StringIO

def commentText(text):
    outL=""
    s = StringIO.StringIO(text)
    for line in s:
        outL += '#  ' + line
    return outL

def procIDtranslate(procid):
    if procid == 0:
        return "total"
    else: 
        return procid

def loadATOMstatInfo(ana_filename):
    try:
        stream = open(ana_filename, 'r')
    except:
        logging.error("Invalid filename: " + ana_filename)
    try:
        AtomStatData = yaml.load(stream,Loader=yaml.CLoader) # Loader=yaml.CLoader
        stream.close()
        return AtomStatData
    except:
        AtomStatData = yaml.load(stream) # Loader=yaml.CLoader
        logging.warning("Using slower python yaml reader. Please install libyaml for faster processing.")
        stream.close()
        return AtomStatData

def getAtomFinalstate(atomYML,procID):
    retFinalstate = ""
    for sbData in atomYML["Sub Processes Information"]:
        if sbData["Process ID"] == procID:
            retFinalstate =  sbData["Process ID Final State"]
    return retFinalstate

def get_results_atom(anaStat, atomOut, xsecin , xsecSubProc):
    results = {}
    warning_list = []
    atomPrint = []
    outData = []

    excluded = ""
    excludedFlags=""

    for atomAnaName in atomOut["Analyses"].keys():
        
        # in pb
        integrated_lumi = atomOut["Analyses"][atomAnaName]["Luminosity"]["Value"] 


        for signal_region in atomOut["Analyses"][atomAnaName]["Efficiencies"]:
            logeffdict={} 
            
            sr_name = signal_region["Efficiency Name"]           
            
            try:
                # print anaStat[atomAnaName]
                dummy =  anaStat[atomAnaName]["SignalRegions"][sr_name] # check if it exists, won't use dummy anymore

                for effData in signal_region["Data"]:
                    excluded = "" 
                    effvalue = effData["Efficiency Value"]
                    efferror = effData["Efficiency Stat Error"]
                            
                    procid = effData["Sub-process ID"]
                    if effvalue > 0:
                        
                        nvis = float(effvalue) * float(xsecSubProc[procid]) * integrated_lumi # number of visible events
                        
                        try:
                            if anaStat[atomAnaName]["SignalRegions"][sr_name]["Type"] != "CutAndCount":
                                print "hello1"
                                logging.warning("Signal region is not Cut&Count " + atomAnaName + "  " + sr_name)
                                raise Exception('Not cut&count')
                            obs95CLnevents = float(anaStat[atomAnaName]["SignalRegions"][sr_name]["UpperLimit95obs"])
                            if obs95CLnevents < 0: #Type: CutAndCount
                                raise Exception('Negative Nevents95')
                        except:
                            try:
                                if anaStat[atomAnaName]["SignalRegions"][sr_name]["Type"] != "CutAndCount":
                                    print "hello1"
                                    logging.warning("Signal region is not Cut&Count " + atomAnaName + "  " + sr_name)
                                    raise Exception('Not cut&count')
                                obs95CLnevents = float(anaStat[atomAnaName]["SignalRegions"][sr_name]["UpperXsecFB95obs"]) * integrated_lumi
                                if obs95CLnevents<0 :
                                    raise Exception('Negative UpperXsecFB95obs')
                            except:
                                logging.warning("Could not find signal region limit information for " + atomAnaName + "  " + sr_name)

                        nvisOn95 = float(nvis) / obs95CLnevents

                        fstate = str(getAtomFinalstate(atomOut,procid))

                        nviserror = [ errorelement  * float(xsecSubProc[procid]) * integrated_lumi for errorelement in efferror ]
                        nvisOn95error = [ errorelement  * float(xsecSubProc[procid]) * integrated_lumi  / obs95CLnevents for errorelement in efferror ]

                        printcolor = bcolors.ENDC

                        if nvisOn95 > 1:
                            excluded = " <- excluded" + bcolors.ENDC
                            excludedFlags = "[Excluded]"
                            printcolor = bcolors.FAIL
                                
                        if nvisOn95 > 1 and (nvisOn95 + nvisOn95error[0]) < 1:  # not enough statistics to claim an exclusion, 1 sigma error is below limit
                            excluded = " excluded? (Low statistics! more MC events" + +bcolors.ENDC
                            excludedFlags = "[LowStatistics,Excluded]"

                        printLine =  [ printcolor + atomAnaName , sr_name,effvalue, nvis, nvisOn95 ]  

                        
                        if errorSHOW:
                            
                            printLine =  [ printcolor + atomAnaName , sr_name, effvalue, efferror[0],efferror[1], nvis, nviserror[0],nviserror[1] , nvisOn95  ]
                            
                        if subprocSHOW:
                            printLine.append( procIDtranslate(procid) )
                            if finalstateSHOW:
                                printLine.append( str(getAtomFinalstate(atomOut,procid)) )
                            printLine.append( excluded )

                            outLine =  [ atomAnaName , sr_name ,effvalue, float(effvalue)*float(xsecSubProc[procid])*  integrated_lumi, float(effvalue)*float(xsecSubProc[procid])*integrated_lumi/float(obs95CLnevents), procIDtranslate(procid), str(getAtomFinalstate(atomOut,procid)), excluded ]
                        
                            atomPrint.append([str(i) for i in printLine])
                            outData.append([str(i) for i in outLine])

                        elif not subprocSHOW and procid == 0:  
                            printLine.append( excluded )
                            outLine =  [ atomAnaName , sr_name ,effvalue, float(effvalue)*float(xsecSubProc[procid])*  integrated_lumi, float(effvalue)*float(xsecSubProc[procid])*integrated_lumi/float(obs95CLnevents), procid,  excluded ]
                            atomPrint.append([str(i) for i in printLine])
        
                            outData.append([str(i) for i in outLine])
        
                        indexEFF =0
                        for effS in atomOut["Analyses"][atomAnaName]["Efficiencies"]:
                            
                            if effS["Efficiency Name"] == sr_name:
        
                                indexPID = 0 
                                for procs in effS["Data"]:                                    
                                    if procs["Sub-process ID"] == procid:
                                        # include it in the ATOM output yaml, we've passed it as a pointer, so the changes affect the input dictionary
                                        atomOut["Analyses"][atomAnaName]["Efficiencies"][indexEFF]["Data"][indexPID]["NvisibleEvents"] = nvis 
                                        atomOut["Analyses"][atomAnaName]["Efficiencies"][indexEFF]["Data"][indexPID]["Nvisible/N95"] = nvisOn95                                        
                                        atomOut["Analyses"][atomAnaName]["Efficiencies"][indexEFF]["Data"][indexPID]["NvisibleEvents Stat Error"] = nviserror
                                        atomOut["Analyses"][atomAnaName]["Efficiencies"][indexEFF]["Data"][indexPID]["Nvisible/N95 Stat Error"] = nvisOn95error
                                        atomOut["Analyses"][atomAnaName]["Efficiencies"][indexEFF]["Data"][indexPID]["Limit Flags"] = excludedFlags                                         
                                        break
                                    indexPID +=1
                            indexEFF += 1
                        excluded = bcolors.ENDC

                
            except KeyError:
                logging.warning("Could not find signal region limit information for " + atomAnaName + "  " + sr_name)              
       
    printHeader = ['Analysis','Signal Region', 'Efficiency' ]

    headErr = ['-','+']

    if errorSHOW:
        printHeader.extend(headErr)

    printHeader.append('Nvis')

    if errorSHOW:
        printHeader.extend(headErr)

    printHeader.extend(['Nvis/N95','Process-ID'])

    if finalstateSHOW:
        printHeader.append('Final State')

    printHeader.append('')     


   # print tabulate(atomPrint, headers="firstrow", tablefmt="rst")
   # atomPrint = sorted(atomPrint, key=operator.itemgetter(4))

    atomPrint.insert(0,printHeader)
    outData.insert(0,printHeader)
    print tabulate(atomPrint, headers="firstrow", tablefmt="rst")
    print bcolors.ENDC
    outData=tabulate(outData, headers="firstrow", tablefmt="rst")
    return outData

# print "\n" + bcolors.OKGREEN + anaNAME + bcolors.ENDC + \
#                     "\n" + bcolors.OKBLUE + "Cuts" +  bcolors.ENDC
#         atomout.printCuts(anaNAME)
#         print "\n"  + bcolors.OKBLUE + "Efficiencies" +  bcolors.ENDC
#         atomout.printEfficiencies(anaNAME)
#         print "\n\n"

import argparse 

if __name__ == "__main__":

    logging.basicConfig(format='%(levelname)s:  %(message)s', level=logging.INFO)

    options = { # 'fastlimdir' : os.path.dirname(os.path.realpath(__file__)), \
               'version_major' : ((__version__).split("."))[0], 
               'version_minor' : ((__version__).split("."))[1], 
               'authors' : __author__,
               'update_interval' : 7, # checking updates every 7 days
               'max_mass' : 2000}

    parser = argparse.ArgumentParser(   # AW added minimal command-line parsers to keep compatibility with scripts
                description='Atom Limit: Automatic Limits Of Models. Run atom results against visible cross section limits.')

    parser.add_argument("-c","--cross-section", dest="tot_cross_section_in", metavar="cross-section in pb",
                       help="Total cross-section", type=float)

    parser.add_argument("-o","-H", "--output", dest="outfilename", default=[],
                       help="write Atom Limit output to file", metavar="BATCHFILE")

    parser.add_argument('-sp','--sub-processes',dest='subprocINFO',action='store_true', 
                     help="Show sub-processes.")
    parser.set_defaults(subprocINFO=False)

    parser.add_argument('-fs','--final-states',dest='finalstateINFO',action='store_true', 
                     help="Show final state information.")
    parser.set_defaults(finalstateINFO=False)

    parser.add_argument('-e','--errors',dest='errorINFO',action='store_true', 
                     help="Show statistical error information.")
    parser.set_defaults(errorINFO=False)

    parser.add_argument("inputfile", 
                      default=[], 
                      help="Atom result file in <yml> format.")

    args = parser.parse_args()
    #print_logo(options)


    analysisInfoDir = "AnaInfo/"

    # atom_share = subprocess.Popen(["atom-config", "--datadir"], stdout=subprocess.PIPE).communicate()[0].strip()
    # atom_version = subprocess.Popen(["atom-config", "--version"], stdout=subprocess.PIPE).communicate()[0].strip()
    # fastlimdir= atom_share+'/Atom-'+atom_version+'/fastlim_dir/'

    
    # Defining a list of the analyses that you want to consider 
    # Information about the analyses is obtained from the files in analyses_info.             
    ana_dict = {}

    # @TODO
    # Fix this to read the ATLAS/CMS infos ... 

    allAnaStat = dict()
    logging.info('Reading ATLAS/CMS limits...') 
    if os.path.exists(analysisInfoDir): 
        ana_list = os.listdir(analysisInfoDir)
        #print ana_list
        for ana_filename in ana_list:
            if ana_filename.endswith(".limit"):
                oneAna=loadATOMstatInfo(os.path.join(analysisInfoDir, ana_filename))
                allAnaStat[oneAna["Name"]] = oneAna
    else:
        ana_list = []
    
    logging.info('Reading Atom file <' + args.inputfile + ">...")    

    try:
        stream = open(args.inputfile, 'r')
    except:
        logging.error("Invalid filename.")
    try:
        AtomData = yaml.load(stream,Loader=yaml.CLoader) # Loader=yaml.CLoader
    except:
        AtomData = yaml.load(stream) # Loader=yaml.CLoader
        logging.warning("Using slower python yaml reader. Please install libyaml for faster processing.")

    atomXsec = 0
    atomXsecSubProc = dict()

    for xsecElement in AtomData["Cross Sections"]:    
        try:
            if xsecElement["Process ID"] == 0:
                atomXsec = xsecElement["Cross Section"]  # Atom gives pb cross-section, we use pb here
                atomXsecSubProc[0] = atomXsec 
            else: 
                atomXsecSubProc[xsecElement["Process ID"]] = xsecElement["Cross Section"]
        except:
            pass

    tot_cross_section_in = 0.
    
    if args.tot_cross_section_in:
        tot_cross_section_in = args.tot_cross_section_in
        xsecInfo = 'Using provided cross-section ' + str(tot_cross_section_in) + ' pb' 
        logging.info(xsecInfo) 
        if args.subprocINFO == True:
            logging.warning("Switching off sub-processes because we do not have xsec information." )
            args.subprocINFO = False 
    elif atomXsec > 0:        
        tot_cross_section_in = atomXsec 
        xsecInfo = 'Taking cross-section from AtomFile: ' + str(tot_cross_section_in) + ' pb ' 
        logging.info(xsecInfo ) 
    else:
        logging.error("Invalid cross-section.")

    AtomRunInfo = yaml.dump(AtomData["Atom Run Info"], default_flow_style=False) 
    
    


    results_8 = {}
    results_7 = {}
   
    ###############################################################


    if args.subprocINFO:
        logging.info('Will show sub-processes')
    if args.finalstateINFO:
        logging.info('Will show final state information')
    if args.errorINFO:
        logging.info('Will show efficiency errors')


    # set flags

    subprocSHOW = args.subprocINFO 
    finalstateSHOW = args.finalstateINFO
    errorSHOW = args.errorINFO

    outData = get_results_atom(allAnaStat, AtomData, tot_cross_section_in,atomXsecSubProc)
       
    # AtomRunInfo = "Input File: "+ args.inputfile + " \n" + xsecInfo + " \n\n" + AtomRunInfo

    # AtomRunInfo= commentText(AtomRunInfo)
    # fout.write(AtomRunInfo)
    # fout.write(outData)
    # fout.close()

   # outputfile="atom_limit.out"
    if args.outfilename:
        outputfile=args.outfilename
        logging.info('Writing result to: <'+ outputfile + '>' )

        yamlName = args.inputfile + ".limit"
        logging.info('Updating '+  args.inputfile  +' including limit result to: <'+ yamlName + '>' )
        with open(yamlName, 'w') as outfile:
            outDump =  yaml.dump(AtomDataInclStat,default_flow_style=False)
            outfile.write( outDump )
            outfile.close()

    logging.info('Writing mathematica input file including limit result to: <'+ args.inputfile + '.mathematica>' )    

    with open(args.inputfile+".mathematica", 'w') as outfile:
         json.dump(AtomData, outfile, indent=4) 

    outfile.close()
    



    exit()

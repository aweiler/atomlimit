#! /usr/bin/env python

"""
"""

__author__ = "A.Weiler"
__version__ = "0.1"


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

def get_results_atom(ana_dict, atomOut, xsecin , xsecSubProc):
    results = {}
    warning_list = []
    atomPrint = []
    outData = []
    print 

    for ana, anaData in ana_dict.iteritems():
        for sr_name, srData in anaData.sr_dict.iteritems():
            logeffdict={}            
            
            atomAnaName = ana
                
            try:
                #print atomAnaName
                for effItem in atomOut["Analyses"][atomAnaName]["Efficiencies"]:
                    #
                    excluded = ""
                    effvalue = 0.
                    #efferror = effItem["Error"][0]["Stat"][1]                   
                    if effItem["Efficiency Name"] == srData.sr_info:
                        #print "found one" , ana, srData.sr_info, effItem["Data"]
                        for effData in effItem["Data"]:
                            effvalue = effData["Efficiency Value"]
                            efferror = effData["Efficiency Stat Error"]
                            #print efferror
                            procid = effData["Sub-process ID"]

                            printLine = []
                            outLine = []

                            if effvalue > 0.:
                                
                                printcolor = bcolors.ENDC
                                

                                nvis = float(effvalue) * float(xsecSubProc[procid]) * anaData.lumi
                                nvisOn95 = float(nvis) / float(srData.UL_nvis_obs)
                                fstate = str(getAtomFinalstate(atomOut,procid))

                                nviserror = [ errorelement  * float(xsecSubProc[procid]) * anaData.lumi for errorelement in efferror ]
                                nvisOn95error = [ errorelement  * float(xsecSubProc[procid]) * anaData.lumi  / float(srData.UL_nvis_obs) for errorelement in efferror ]

                                #print srData.UL_nvis_obs, nvisOn95, nvisOn95error, (nvisOn95 + nvisOn95error[0])

                                if nvisOn95 > 1:
                                    excluded = " <- excluded"
                                    printcolor = bcolors.FAIL
                                
                                if nvisOn95 > 1 and (nvisOn95 + nvisOn95error[0]) < 1: # 
                                    excluded = " <- excluded? (LOW STAT)"
                                    #print excluded
                                    printcolor = bcolors.WARNING

                                printLine =  [ printcolor + ana , srData.sr_info,effvalue, nvis, nvisOn95,  ]

                                if errorSHOW:
                                    #printLine =  [ printcolor + ana , srData.sr_info, effvalue, str(efferror), nvis, str(nviserror) , nvisOn95  ]
                                    printLine =  [ printcolor + ana , srData.sr_info, effvalue, efferror[0],efferror[1], nvis, nviserror[0],nviserror[1] , nvisOn95  ]
                                    #print printLine

                                if subprocSHOW:

                                    printLine.append( procid )
                                    if finalstateSHOW:
                                        printLine.append( str(getAtomFinalstate(atomOut,procid)) )

                                    printLine.append( excluded )

                                    outLine =  [ ana , srData.sr_info,effvalue, float(effvalue)*float(xsecSubProc[procid])*  anaData.lumi, float(effvalue)*float(xsecSubProc[procid])*anaData.lumi/float(srData.    UL_nvis_obs), procid, str(getAtomFinalstate(atomOut,procid)), excluded ]
                                    atomPrint.append([str(i) for i in printLine])
                                    outData.append([str(i) for i in outLine])    

                                elif not subprocSHOW and procid == 0:    

                                    printLine.append( excluded )
                                    outLine = [ ana , srData.sr_info,effvalue, float(effvalue)*float(xsecin)*   anaData.lumi, float(effvalue)*float(xsecin)*anaData.lumi/float(srData. UL_nvis_obs), procid, excluded ]
                                    atomPrint.append([str(i) for i in printLine])
                                    outData.append([str(i) for i in outLine])
    
                                # elif subprocSHOW and not finalstateSHOW:

                                #     printLine.append( procid )
                                #     printLine.append( excluded )

                                #     outLine =  [ ana , srData.sr_info,effvalue, float(effvalue)*float(xsecSubProc[procid])*  anaData.lumi, float(effvalue)*float(xsecSubProc[procid])*anaData.lumi/float(srData.    UL_nvis_obs), procid, excluded ]
                                #     atomPrint.append([str(i) for i in printLine])
                                #     outData.append([str(i) for i in outLine])
                                #     #atomPrint.append(printLine)
                                
                                # print "debug opts ",subprocSHOW , procid , finalstateSHOW

                                
                                excluded = ""
            except KeyError:
                #print "Unexpected error:", sys.exc_info()[0]
                pass    

            except:
                print "Unexpected error:", sys.exc_info()
                pass    
            #print "done"       

    # printHeader = ['Analysis','Signal Region', 'efficiency', 'Nvis','Nvis/N95','Process-ID','']

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

    atomPrint.insert(0,printHeader)
    outData.insert(0,printHeader)
    print atom_limit.tabulate.tabulate(atomPrint, headers="firstrow", tablefmt="rst")
    print bcolors.ENDC
    outData=atom_limit.tabulate.tabulate(outData, headers="firstrow", tablefmt="rst")
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

    parser.add_argument("-c","--cross-section", dest="tot_cross_section_in", metavar="cross-section in fb",
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
        print ana_list
        for ana_filename in ana_list:
            oneAna=loadATOMstatInfo(os.path.join(analysisInfoDir, ana_filename))
            allAnaStat[oneAna["Name"]] = oneAna
    else:
        ana_list = []

    # if os.path.exists(os.path.join(fastlimdir,'7TeV')):         
    #     ana_list_7 = os.listdir(os.path.join(fastlimdir,'7TeV'))
    #     ana_dict_7 = atom_limit.read_data.get_ana_dict(fastlimdir,ana_list_7, 7)
    #     ana_dict.update(ana_dict_7)
    # else:
    #     ana_list_7 = []
    print allAnaStat
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
                atomXsec = xsecElement["Cross Section"] * 1000 # Atom gives pb cross-section, we use fb
                atomXsecSubProc[0] = atomXsec * 1000
            else: 
                atomXsecSubProc[xsecElement["Process ID"]] = xsecElement["Cross Section"] *1000 
        except:
            pass

    tot_cross_section_in = 0.
    
    if args.tot_cross_section_in:
        tot_cross_section_in = args.tot_cross_section_in
        xsecInfo = 'Using provided cross-section ' + str(tot_cross_section_in) + ' fb' 
        logging.info(xsecInfo) 
        if args.subprocINFO == True:
            logging.warning("Switching off sub-processes because we do not have xsec information." )
            args.subprocINFO = False 
    elif atomXsec > 0:        
        tot_cross_section_in = atomXsec 
        xsecInfo = 'Taking cross-section from AtomFile: ' + str(tot_cross_section_in) + ' fb ' 
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

    outData = get_results_atom(ana_dict_8, AtomData, tot_cross_section_in,atomXsecSubProc)
    
    outputfile="atom_limit.out"
    if args.outfilename:
        outputfile=args.outfilename
    logging.info('Writing result to: <'+ outputfile + '>' )

    fout = open(outputfile, "w")
    AtomRunInfo = "Input File: "+ args.inputfile + " \n" + xsecInfo + " \n\n" + AtomRunInfo

    AtomRunInfo= commentText(AtomRunInfo)
    fout.write(AtomRunInfo)
    fout.write(outData)
    fout.close()
    exit()
import os
import re
import shutil


selection = open("selection.cpp","r")
selectionLines = selection.readlines()
selection.close()

settings = open("settings.cpp","r")
settingsLines = settings.readlines()
settings.close()

spectrumPosition=["all","onPeak","offPeak","farOffPeak1","farOffPeak2"]
type=["signal","background","sigBg"]

for cutType in spectrumPosition:
    for theType in type:
        print "will do", cutType, theType
        os.mkdir("ggDir_"+cutType+"_"+theType)
        shutil.copy("selection.h","ggDir_"+cutType+"_"+theType+"/selection.h")
        outFile = open("ggDir_"+cutType+"_"+theType+"/selection.cpp","w")
        for line in selectionLines:
            outFile.write(line)
            if "isMassCut" in line:
                if (cutType=="onPeak"):
                    outFile.write("require(llll().m() > 110.);\n")
                    outFile.write("require(llll().m() < 140.);")
                elif (cutType=="offPeak"):
                    outFile.write("require(llll().m() > 140.);")
                elif (cutType=="farOffPeak1"):
                    outFile.write("require(llll().m() > 220.);\n")
                    outFile.write("require(llll().m() < 300.);")
                elif (cutType=="farOffPeak2"):
                    outFile.write("require(llll().m() > 300.);")
        outFile = open("ggDir_"+cutType+"_"+theType+"/settings.cpp","w")
        for line in settingsLines:
            if ("CONTINUUMBACKGROUND" in line) and (theType=="signal"):
                continue
            if ("HIGGSSIGNAL" in line) and (theType=="background"):
                continue
            outFile.write(line)


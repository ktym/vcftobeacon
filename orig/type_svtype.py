#!/usr/bin/python

import sys, getopt
import pandas as pd
import sys

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile="])
   except getopt.GetoptError:
      print ('type_svtype.py -i <inputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('type_svtype.py -i <inputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      print ('Input file is "', inputfile)
      data = pd.read_csv(inputfile,sep = "\t", skiprows=[0], header=None, names=["TYPE", "SVTYPE"], na_values = ".")
      #if the number of NAs is 0 for both columns, SVTYPE column will be deleted (case 2):
      if data["SVTYPE"].isnull().sum() == data["TYPE"].isnull().sum() == 0:
         data=data.drop(columns="SVTYPE")
         data["TYPE"].to_csv("type", sep="\t", header=True, index = False, na_rep = ".")
      #elif TYPE and SVTYPE is all of NAs, we save the file with ".":
      elif data["TYPE"].isnull().sum() == data["SVTYPE"].isnull().sum() == len(data):
         data["TYPE"].to_csv("type", sep="\t", header=True, index = False, na_rep = ".")
      #elif SVTYPE only is all NAs we save the column of TYPE:
      elif data["SVTYPE"].isnull().sum() == len(data):
         data["TYPE"].to_csv("type", sep="\t", header=True, index = False, na_rep = ".")
      #elif TYPE only is all NAs we save the column of SVTYPE:
      elif data["TYPE"].isnull().sum() == len(data):
         data=data.drop([TYPE])
         data.rename(columns={'SVTYPE':'TYPE'}, inplace=True)
         data["TYPE"].to_csv("type", sep="\t", header=True, index = False, na_rep = ".")

if __name__ == "__main__":
   main(sys.argv[1:])
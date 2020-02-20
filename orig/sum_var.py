#!/usr/bin/python

import sys, getopt
import pandas as pd

def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print ('sum_var.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('sum_var.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      print ('Input file is "', inputfile)
      data = pd.read_csv(inputfile,sep = "\t", header=0, na_values = ".")
      data["SUM"]=data["nHet"] +data["nHomAlt"]
      data["SUM"].to_csv("Sum_variants", sep="\t", header=True, index = False)


if __name__ == "__main__":
   main(sys.argv[1:])
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 13 12:49:30 2018

@author: sean.montgomery
"""
import sys

def read_file(filename,binsize):
  """(string) -> string
  
  Opens filename in the current directory.
  """
  
  with open(filename) as file:
      output = "" #Sets an empty list
      binsize = int(binsize)
      binstep = binsize
      bin = binsize
      binsum = 0
      bincount = 0
      for line in file: #Cycles through lines in the file
          line=line.split()
          if line[0] == "chr1":
              chr = line[0]
          if line[0] == chr:
              if int(line[1]) < bin:
                  binsum += float(line[3])
                  bincount += 1
              elif int(line[1]) >= bin:
                  binavg = binsum/bincount
                  output += chr + '\t' + str(int(bin-binstep/2)) + '\t' + str(binavg) + '\n'
                  bin += binstep
                  binsum = float(line[3])
                  bincount = 1
          else:
              chr = line[0]
              binsum = float(line[3])
              bincount = 1
              bin = binsize
      #output += chr + '\t' + str((bin-binstep)+(binstep-(bin-int(line[1])))/2) + '\t' + str(binsum/bincount) + '\n'
  return output
  
def write_file(output, new_filename):           
    file = open(new_filename, 'w') #Creates a new file with the name $new_filename
    file.write(output) #Writes its input into the new file
    file.close() #Closes the new file
    
    
if __name__ == '__main__':
    file_name = sys.argv[1]
    binsize = sys.argv[2]
    newfile = sys.argv[3]
    new_bed=read_file(file_name,binsize)
    write_file(new_bed,newfile)
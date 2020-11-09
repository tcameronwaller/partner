# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 11:26:23 2016

@author: Cameron
"""

#Author: T Cameron Waller
#First edit: 22 February 2016
#Last edit: 14 June 2016

#Description:
#QuanLynx exports data in a text file using the "Complete Summary" Export option.
#This program parses the text file and exports the data to a tab-delimited table

###########################################################################
# Definitions of functions

def read(list_row=None, file=None):
    
    # Documentation string
    
    """
    This function reads information from QuanLynx Complete Summary.
    """
    
    # Procedure
    
    if list_row is not None and file is not None:

        with open(file, "r") as file_in:
            
            #Print file lines
            #print(file.read())
            #print(file.readline())
            #print(file.readlines())
            
            #Iterate over lines in file
            for line in file_in:
                
                #print(iteration)
                
                #Identify and store analyte from block header in file
                if line.strip()[0:8] == "Compound":
                    #print(line)
                    list_line_split = line.strip().split()
                    analyte = list_line_split[2].strip()
                    #print(analyte)
                    
                # Identify and store column header names from first header
                # line in file
                elif not list_row and line.strip()[0:1] == "#":
                    #print(line)
                    list_line_split = line.strip().split("\t")
                    #print(list_line_split[1])
                    #print(list_line_split[1:11])
                    list_column = ["Analyte"]
                    list_column.append(list_line_split[1])
                    list_column.append(list_line_split[4])
                    list_column.append(list_line_split[5])
                    #print(list_column)
                    #list_row.append(list_column[:])
                    delimiter = "\t"
                    list_row.append(delimiter.join(list_column[:]))
                    
                #Identify and store values from lines in file
                elif line.strip().split("\t")[0].isdigit():
                    #print(line)
                    list_line_split = line.split("\t")
                    #print(list_line_split[2])
                    #print(list_line_split[5])
                    #print(list_line_split[6])
                    #print(list_line_split[2:12])
                    list_column = [analyte]
                    list_column.append(list_line_split[2])
                    list_column.append(list_line_split[5])
                    list_column.append(list_line_split[6])
                    #print(list_column)
                    #list_row.append(list_column[:])
                    delimiter = "\t"
                    list_row.append(delimiter.join(list_column[:]))
                #iteration = iteration + 1
                
        return list_row

#Print list of rows
#print(list_row)
#print(list_row[0])
#print(list_row[1])
#print(list_row[2])
#print(list_row[3])
#print(list_row[4])
#print(list_row[5])


###########################################################################
# Script

#Note:

#I should eventually put this code in a function that accepts directory, in file, and out file as arguments
#It would also be nice if the master/main script automatically parsed files
# for both standard and apex
#Maybe I should make the function able to create an empty out file itself so that I don't have to do that separately
#Then I should also put this function in a module or something

#Specify directory and in and out files
import os
#directory = os.path.join(os.path.expanduser("~"), "6521", "test", "cameron")
directory = os.path.join(
    "C:\\", "Data", "Local", "Research_Rutter", "Projects", "2016_Schell",
    "Mass_Spectrometry_2016-08-15", "Quantity")
print(directory)
#directory = "C:\Data\Research_Rutter\Projects\Metabolite_Analysis\\"
#file_in_name = os.path.join(directory,
#                            "Quantity_Standard_Complete_Summary.txt")

file_in_name_1 = os.path.join(directory,
                            "quantity_standard_1_complete_summary.txt")
print(file_in_name_1)
file_in_name_2 = os.path.join(directory,
                            "quantity_standard_2_complete_summary.txt")
print(file_in_name_2)
file_out_name = os.path.join(directory, "quantity_standard_report.txt")

#####

file_in_name_1 = os.path.join(directory,
                            "quantity_apex_1_complete_summary.txt")
print(file_in_name_1)
file_in_name_2 = os.path.join(directory,
                            "quantity_apex_2_complete_summary.txt")
print(file_in_name_2)
file_out_name = os.path.join(directory, "quantity_apex_report.txt")
print(file_out_name)

#Read information from file
#Note, as written, the parser stores individual rows as delimited strings
#This is to simplify writing the information to file
#If I wanted to use the information in the program, then I would instead store the row information as a list
#The total data set would be stored in a list of lists
#The list of rows refers to lists of columns
#iteration = 1
list_row = []
list_row = read(list_row=list_row, file=file_in_name_1)
list_row = read(list_row=list_row, file=file_in_name_2)

#Write information to file

with open(file_out_name, "w") as file_out:
    
    #Iterate over lines in file
    for row in list_row:
        
        #Write row as line in file
        #print(row)
        file_out.write(row + "\n")
#print("Is file_out closed?", file_out.closed)








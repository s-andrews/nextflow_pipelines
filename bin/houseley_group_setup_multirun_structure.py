#!/usr/bin/env python3

import os
import re 
import sys
import argparse
import csv

# Script to help organise the processing of run folders with multiple sequencing types within them
    # First argument is the run folder name
    # Second argument is the run_info.csv 
        # this is now a set format for the Housley group
            # Sample name, Barcode 1, Barcode 2 (if using), unique library group identifier, Type of run, Genome, Trael indexes (if using), Trael names (if using)
    # optional argument is --dual
        # default behaviour is for single index, --dual specifies dual index
    # output directories/symlinks are created in current working directory
            
# More detail on columns
    # Sample name, we ignore this column, just for ease of use for the group
    # Barcode 1/2 are needed to identify the samples
    # unique library group identifier, this is a identifier for groups of samples so that they can easily be identified in case of re-processing
    # Type of run, what kind of nf pipeline to use 
    # Genome, what kind of genome to use
    # Trael indexes, a comma separated list of indexes to include, if no-index needed this will be stated as "no index"
    # Trael names, a comma separated list of names for each index, must be in the same order as indexes. If "no index" leave blank


# The run descriptor columns (everything apart from barcode/sample name) will be joined to create an informative directory name
    # These directories are created within the current working directory
    # Within each of these sub-directories
        # Fastq files whose barcodes match the run descriptions are stored as symlinks
        # Fastq file symlink names include the unique library identifier (added after existing sample name)

def main():

    # read in options from the command line
    options = parse_arguments()

    run_folder = options.run_folder
    info_file = options.run_info_file
    no_header= options.no_header

    if not options.run_info_file.endswith(".csv"):
        parser.error("The run info file must be a CSV file.")

    #set working directory (place directory structure will be made and symlinks will be stored)
    working_dir = os.getcwd()
    
    #set the path to the seqfac fastq files
    seqfac_path="/bi/seqfac/seqfac/"+run_folder+"/Unaligned/Project_External/Sample_lane1/"
    
    # create the directory structure in the current working directory
    run_dict = read_run_info(info_file,no_header)

    # create links to fastq files within the appropriate sub-directories 
    make_multi_run_structure(run_dict,working_dir,seqfac_path)

def make_multi_run_structure(run_dict, working_dir,seqfac_path):
    for run_type, barcodes in run_dict.items():
        dir_path = os.path.join(working_dir, run_type)
        make_dir(dir_path)
        make_links(seqfac_path, dir_path, run_type,barcodes)

def make_links(seqfac_path, dest_dir, run_type,barcodes):
    
    #the barcodes to search for
    # added a fix here to handle lower case bases in barcode sheet
    patterns = [re.compile(barcode.upper()) for barcode in barcodes]

    # Extract the unique ID from the dir_name (i.e. run_type)
    ID = run_type.split('_')[0]
    
    # for the fastq files in the seqfac run folder look for a barcode match and create link in subdirectory
    for root,dirs,files in os.walk(seqfac_path):
        for file in files:
            if not file.endswith(".fastq.gz"):
                continue
            file_path = os.path.join(root, file)
            
            if any(pattern.search(file) for pattern in patterns):
                ending = r"(_L\d{3}_R\d\.fastq\.gz)$"
                # create the link name from the first entry in barcodes which is the reduced run_info
                link_name = re.sub(ending, f"_{ID}\\1", file)
                link_path = os.path.join(dest_dir,link_name)

                if not os.path.exists(link_path):
                    os.symlink(file_path, link_path)
                else:
                    print(f"The link already exists: {link_path}")


def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def read_run_info(file,no_header):
    
    run_dict = {}
    

    with open(file, mode='r', newline='') as infile:
        reader = csv.reader(infile)
        #they dont always put a header line
        if not no_header:
            next(reader)    
        
        for row in reader:
            #remove empty cells & then also check for commas or spaces within cells
            line = [format_cell(cell) for cell in row if format_cell(cell)]
            
            if check_if_dual(line[2]):
                dual_barcode = True
            else:
                dual_barcode = False

            # check if the barcode is dual indexed if so join these together
            if dual_barcode:
                dir_name = '_'.join(line[3:])
                barcode = '_'.join(line[1:3])
            else:
                dir_name = '_'.join(line[2:])
                barcode = line[1]
                

            if dir_name in run_dict:
                run_dict[dir_name].append(barcode)
            else:
                run_dict[dir_name] = [barcode]

    return(run_dict)

def check_if_dual(second_entry):
    pattern = r'^[CATG]*$'
    
    # Check if the string matches the pattern
    if re.fullmatch(pattern, second_entry):
        return True
    else:
        return False

def format_cell(cell):
# function to format the contents of a cell to fix any issues that might mess with naming
    cell = cell.strip()

    if ',' in cell:
        items = cell.split(',')
        items = [item.strip() for item in items]  # Remove leading/trailing spaces
        cell = '-'.join(items)

    #add a check to remove spaces from cells
    cell = cell.replace(" ", "_")

    #add a check to remove brackets from cells
    cell = cell.replace('(', '').replace(')', '')
    return cell

def parse_arguments():

    parser = argparse.ArgumentParser(description="Create a Directory Structure for Multi-sequencing runs & link to fastq files")
    #parser.add_argument("--dual", action='store_true', default=False)
    parser.add_argument("run_folder", help="Name of the run folder")
    parser.add_argument("run_info_file", help="Path to the run info CSV file")
    parser.add_argument("--no_header", action='store_true', default=False, help="specify if there is no header in file, default is false")

    return parser.parse_args()
      
if __name__ == "__main__":
    main()
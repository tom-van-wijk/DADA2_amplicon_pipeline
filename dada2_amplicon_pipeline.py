#!/usr/bin/env python


# Name:		dada2_amplicon_pipeline.py
# Author:	Tom van Wijk
# Date:		21-06-2018
# Licence:	GNU General Public License v3.0 (copy provided in directory)

# For detailed information and instruction on how to install and use this software
# please view the included "README.md" file in this repository


# import python libraries
from argparse import ArgumentParser
from time import gmtime, strftime
import os
import sys
import logging
import logging.handlers
import random


# Function to parse the command-line arguments
# Returns a namespace with argument keys and values
def parse_arguments(args, log):
	log.info("Parsing command line arguments...")
	parser = ArgumentParser(prog="dada2_amplicon_pipeline.py")
	parser.add_argument("-i", "--indir", dest = "input_dir",
		action = "store", default = None, type = str,
		help = "Location of input directory (required)",
		required = True)
	parser.add_argument("-o", "--outdir", dest = "output_dir",
		action = "store", default = "inputdir", type = str,
		help = "Location of output directory (default=inputdir)")	
	parser.add_argument("-x", "--savetemp", dest = "save_temp",
		action = "store", default = "false", type = str,
		help = "Option to save temporary files (default=false)")
	return parser.parse_args()


# Function creates logger with handlers for both logfile and console output
# Returns logger
def create_logger(logid):
	# create logger
	log = logging.getLogger()
	log.setLevel(logging.INFO)
	# create file handler
	fh = logging.FileHandler(str(logid)+'_dada2_amplicon_pipeline.log')
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(logging.Formatter('%(message)s'))
	log.addHandler(fh)
	# create console handler
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(logging.Formatter('%(message)s'))
	log.addHandler(ch)
	return log
	

# Function creates a list of files or directories in <inputdir>
# on the specified directory depth
def list_directory(input_dir, obj_type, depth):
	dir_depth = 1
	for root, dirs, files in os.walk(input_dir):
		if dir_depth == depth:
			if obj_type ==  'files':
				return files
			elif obj_type == 'dirs':
				return dirs
		dir_depth += 1
		

# Function to run the dada2classification script
# input: directory with processed and filtered/sorted .fastq filepairs
def run_dada2_classify(input_directory, reference_directory):
	os.system("dada2_classify_contigs.R "+input_directory+" "+reference_directory
		+" | tee "+input_directory+"/dada2_classify_contigs.log")


# Function closes logger handlers
def close_logger(log):
	for handler in log.handlers:
		handler.close()
		log.removeFilter(handler)


# MAIN function
def main():
	# create logger
	logid = random.randint(99999, 9999999)
	log = create_logger(logid)
	# parse command line arguments
	args = parse_arguments(sys.argv, log)
	# creating output directory
	if args.output_dir == 'inputdir':
		output_dir = os.path.abspath(args.input_dir+"/dada2_amplicon_pipeline_output")
	else:
		output_dir = os.path.abspath(args.output_dir)
	log.info("Creating output directory: "+output_dir)
	os.system("mkdir -p "+output_dir)
	# create directories for logfiles, temporary data and output in output_dir
	os.system("mkdir -p "+output_dir+"/temp")
	os.system("mkdir -p "+output_dir+"/logfiles")
	os.system("mkdir -p "+output_dir+"/quality_control")
	# copy data to temp dir
	os.system("cp "+args.input_dir+"/*.fastq "+output_dir+"/temp/.")
	os.system("cp "+args.input_dir+"/*.fastq.gz "+output_dir+"/temp/.")
	# run dada2 script
	log.info("Executing the dada2 classification script")
	run_dada2_classify(os.path.abspath(output_dir+"/temp"), os.environ['DADA2_AMPLICON_HOME'])
	# Move relevant dada2 output to output and output/logfiles directory
	os.system("mv "+output_dir+"/temp/*.log "+output_dir+"/logfiles/")
	os.system("mv "+output_dir+"/temp/*.pdf "+output_dir+"/quality_control/")
	os.system("mv "+output_dir+"/temp/*.rds "+output_dir+"/")
	os.system("mv "+output_dir+"/temp/*.txt "+output_dir+"/")
	# Remove the temp directory if temp parameter is not set to anything but false
	if args.save_temp == "false":
                log.info("save_temp parameter is false (default). Removing temporary files and directories...")
		os.system("rm -rf "+output_dir+"/temp")
	else:
		log.info("save_temp parameter is not false (set in input parameter). Keeping temporary files and directories...")
	# close logger handlers
	log.info("Closing logger and finalising dada2_amplicon_pipeline.py...")
	close_logger(log)
	# move logfile to output directory
	os.system("mv "+str(logid)+"_dada2_amplicon_pipeline.log "+output_dir+"/logfiles/dada2_amplicon_pipeline.log")


main()

# -*- coding: utf-8 -*-
# @Author: tamascogustavo
# @Date:   2021-09-27 10:44:15
# @Last Modified by:   gustavotamascohotmail.com
# @Last Modified time: 2021-10-08 10:27:56


#Imports
import pandas
from time import time
import cobra.test
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import os
from sys import argv


#Functions
def help_menu():
	print("""
		usage: python3 simulating_knockouts.py [-h] [--sg | --dg | --sr] [--dr]
             [-tdr] [--cg] [-cr] 

		Perform gene or reaction deletion in a Metabolic Model created by Chimera or similar tool

		positional arguments:
		  INPUT                 Input (txt file with reactions or genes to be deleted).
		                        When used with -gf for gene file.
		                        When used with --rf for reaction file.

		optional arguments:
		  -h, --help            		show this help message and exit
		  --sg <faa_file target.txt>    To perform single gene deletion
		  --dg <faa_file target.txt>   	To perform double gene deletion
		  --sr reaction file        	To perform single reaction deletion
		  --dr reaction file        	To perform double reaction deletion
		  --tdr       					To perform duble reaction deletion in only two specific
		  								reactions
		  --cg                 			To perform single gene deletion for all genes in the model
		  --cr                 			To perform single reaction deletion all reactions in the model
		""")

def control(argument, model):
	try:
		if argv[1] == "--help" or argv[1] == "-h":
			help_menu()
	except:
		pass
	try:	
		if argv[1] =="--sg":
			with open(argv[2], "r") as faa_file:
				genes_info = parse_faa_info(faa_file)
	
			with open(argv[3], "r") as targets_file:
				targets = get_targets(targets_file)
	
			knockout_ids = mine_model_gene_id(targets, genes_info)
			single_gene_deletion_for_all_list(model,knockout_ids)
	except:
		pass
	try:
		if argv[1] =="--dg":
			with open(argv[2], "r") as faa_file:
				genes_info = parse_faa_info(faa_file)
	
			with open(argv[3], "r") as targets_file:
				targets = get_targets(targets_file)
	
			knockout_ids = mine_model_gene_id(targets, genes_info)
			double_gene_deletion_for_all_list(model, knockout_ids)
	except:
		pass
	
	try:
		if argv[1] =="--sr":
			reactions_list =argv[2]
			rec_file = open(reactions_list, "r")
			rec_list = rec_file.readlines()
			rec_list = [x.strip("\n") for x in rec_list]
			single_target_reaction_deletion(model, rec_list)
	except:
		pass
	try:	
		if argv[1] =="--dr":
			reactions_list =argv[2]
			rec_file = open(reactions_list, "r")
			rec_list = rec_file.readlines()
			rec_list = [x.strip("\n") for x in rec_list]
			double_target_reaction_deletion_for_all(model, rec_list)
	except:
		pass
	try:
		if argv[1] =="--tdr":
			reactions_list =argv[2]
			rec_file = open(reactions_list, "r")
			rec_list = rec_file.readlines()
			rec_list = [x.strip("\n") for x in rec_list]
			specific_double_target_reaction_deletion(model, rec_list)
	except:
		pass
	try:
		if argv[1] =="--cg":
			complete_genes_deletion(model)
	except:
		pass
	try:
		if argv[1] =="--cr":
			complete_reactions_deletion(model)
	except:
		pass
def import_model(model_name):
	'''
	This function import a xml or sbml model to cobrapy

	:param model_name: is a string with the name of the file
	:return model: cobrapy model object
	'''
	model = cobra.io.read_sbml_model(model_name, f_replace={})

	return model

def parse_faa_info(faa_file):
	'''
	This function parses a faa file and produces a dict 

	The dict contain the name of the gene and its code

	:param faa_file: is a faa file used to build the model in previous step
	:return a dict: k = is the gene code in the faa file v= is the gene id
	'''
	faa_info = {}
	for line in faa_file:
		if line.startswith(">"):
			if "[gene=" in line:
				index = line.strip().split()[0].split("|")[0].strip(">")
				gene_locus = line.strip().split()[0].split("|")[1]
				gene_locus = "{}_{}".format(index,gene_locus.replace(".","_"))
				gene_id = line.strip().split()[1].split("=")[1][0:-1]
				faa_info[gene_locus] = gene_id
	return faa_info

def get_targets(file):
	'''
	This function parses a txt file and return all the genes in it

	:param file: txt file containing the gene name, one by row
	:ruturn a list of all genes in the file

	'''
	all_targets = []
	for line in file:
		target = line.strip()
		all_targets.append(target)
	return all_targets

def mine_model_gene_id(targets, genes_info):
	'''
	This function is used to merge info from the targets and faa file

	:param targets: list of targets (genes)
	:param genes_info: a dict. k = is the gene code in the faa file v= is the gene id
	:return a new dict k= are the code in faa file and v is the gene id of the target
	'''
	knockout_dict = {}
	for target in targets:
		if target in genes_info.values():
			knock_gene = list(genes_info.keys())[list(genes_info.values()).index(target)]
			knockout_dict[knock_gene] = target
	return knockout_dict

def single_gene_deletion_for_all_list(model,knockout_ids):
	'''
	This function controls the next one
	:param model: is the cobra model
	:param knockout_ids: dict containing info of the targets and its code in the model
	'''
	for target in knockout_ids.keys():
		knock_out_list(model, target, knockout_ids)


def knock_out_list(model, gene, all_info):
	'''
	This function is reponsable to perform knockout in series for all the targets in the dict

	:param model: is the cobrapy model
	:gene: is a str of the gene code in the model
	:all_info: is a dict that is used to translate the results
	'''
	format_gene = "G_{}.knock_out()".format(gene)
	cmd = "model.genes.{}".format(format_gene)
	try:
		print('<<<Complete model>>>: ', model.optimize())
		with model:
		   exec(cmd)
		   print(f'{all_info[gene]} knocked out: ', model.optimize())
	except: 
		print("Gene {} was not incoporated to the model automatically".format(all_info[gene]))
	
	
def double_gene_deletion_for_all_list(model, knockout_ids):
	'''
	This function does double deletion for all combination of genes provided

	:param model: is the cobrapy model
	:param knocout_ids: is a dict with genes ids and their code in the genome
	:return a file with the results
		
	'''
	outname = 'double_gene_deletion_for_all_list.csv'
	if os.path.exists(outname):
		print(f"{outname} already generated")
	else:
		list_gene_objects = []
		for gene in knockout_ids.keys():
			try:
				format_gene = "G_{}".format(gene)
				gene_obj = eval(f"model.genes.{format_gene}")
				list_gene_objects.append(gene_obj)
			except:
				print(f"{gene} not detected in model")
		

		result = (double_gene_deletion(model, list_gene_objects, list_gene_objects, processes = 12).round(4))
		print(f"{outname} was created")
		result.to_csv(outname)



def single_target_reaction_deletion(model, reaction_list):
	'''
	This function does a single deletion of reaction
	Delete all reactions one by one in the txt file

	:param model: is the cobrapy model
	:param reaction_list: is a list with all reactions provided in the txt file
	'''
	for reaction in reaction_list:
		format_reaction = "R_{}.knock_out()".format(reaction)
		cmd = "model.reactions.{}".format(format_reaction)
		try:
			print('<<<Complete model>>>: ', model.optimize())
			with model:
				exec(cmd)
				print(f'{reaction} knocked out: ', model.optimize())
		except: 
			print("Reaction {} was not incoporated to the model automatically".format(reaction))


def specific_double_target_reaction_deletion(model, reaction_list):
	'''
	This function does a double deletion of reactions
	Need to provide a txt file with 2 lines one reaction in each

	:param model: is the cobrapy model
	:param reaction_list: is a list with all reactions provided in the txt file
	'''

	f_reac1 = "R_{}".format(reaction_list[0])
	f_reac2 = "R_{}".format(reaction_list[1])

	reac1_r = eval("model.reactions.{}".format(f_reac1))
	reac2_r = eval("model.reactions.{}".format(f_reac2))
	print(double_reaction_deletion(model, [reac1_r], [reac2_r], processes = 12).round(4))


def double_target_reaction_deletion_for_all(model, reaction_list):
	'''
	This function does a double deletion of reaction 
	Delete all reactions possible double combinations

	:param model: is the cobrapy model
	:param reaction_list: is a list with all reactions provided in the txt file

	:return a csv file with the results
	'''
	outname = 'double_reaction_deletion_for_all_list.csv'
	if os.path.exists(outname):
		print(f"{outname} already generated")
	else:
		list_reac_objects = []
		for reac in reaction_list:
			try:
				format_reac = "R_{}".format(reac)
				reac_obj = eval(f"model.reactions.{format_reac}")
				list_reac_objects.append(reac_obj)
			except:
				print(f"{reac} not detected in model")
		

		result = (double_reaction_deletion(model, list_reac_objects, list_reac_objects, processes = 12).round(4))
		result.to_csv(outname)
		print(f"{outname} was created")
		


def complete_genes_deletion(model):
	'''
	This function do a single deletion for all genes in the model

	:param model: is the cobrapy model
	:return a csv file with the results of each single deletion
	'''
	if os.path.exists("all_single_gene_knockout.csv"):
		print("all_single_gene_knockout.csv already exists")
	else:
		result = (single_gene_deletion(model, model.genes[0:], processes = 12).round(4))
		result.to_csv('all_single_gene_knockout.csv')
		print("all_single_gene_knockout.csv was created")

def complete_reactions_deletion(model):
	'''
	This function do a single deletion for all reactions in the model

	:param model: is the cobrapy model
	:return a csv file with the results of each single deletion
	'''
	if os.path.exists("all_single_reactions_knockout.csv"):
		print("all_single_reactions_knockout.csv already exists")
	else:
		result = (single_reaction_deletion(model, model.reactions[0:], processes = 12).round(4))#
		result.to_csv('all_single_reactions_knockout.csv')
		print("all_single_reactions_knockout.csv was created")

def main():

	try:
		input_data = argv[1]
		if input_data == "-h" or input_data =="--help":
			control(input_data, "help only")
		else:
			model_name = str(input("Enter de model file name that you want to use: "))
			model = import_model(model_name)
			control(input_data, model)
	except:
		print("No input, use -h for help")

	

if __name__ == '__main__':
	main()
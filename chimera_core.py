# -*- coding: utf-8 -*-
# @Author: gustavotamascohotmail.com
# @Date:   2021-10-04 14:05:50
# @Last Modified by:   tamascogustavo
# @Last Modified time: 2021-10-20 10:25:17



# import statements

from sys import argv
import os.path
import subprocess
import sys
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import gapfill
import os
import cobra.manipulation
import cplex
from os import listdir
from os.path import isfile, join
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
from collections import Counter
import escher
from escher import Builder



# functions and classes
def help_menu():
    print("""
        genome_path, model_name, universe, medium, initiate)
        usage: python3 chimera_core.py <organism.faa> <biomass function> <growth medium for gapfill> 
        <medium to initialize model>

        e.g: python3 chimera_core.py e_coli_test.faa gramneg M9 M9
        [-h]

        Perform model construction, evaluation and visualization with Chimera 

        positional arguments:
          INPUT                 Inputs.
                                faa file: e_coli.faa
                                biomass function: gramneg, grampos, bacteria
                                growth medium for gapfill: LB, LB[-O2], M9, M9[-O2], M9[glyc]
                                    obs: More than one can be use. Use "," to separate them with no spaces 
                                medium to initialize: One of the selected medium
                                    obs: will be use for the simulations

        optional arguments:
          -h, --help                    show this help message and exit
        """)

def control(argument):
    try:
        if argv[1] == "--help" or argv[1] == "-h":
            help_menu()
        else:
            '''
            Setting the variable
            '''
            dirpath = os.getcwd()
            escher.rc['never_ask_before_quit'] = True
            all_files = list_files(dirpath)
            genome_file = argument
            genome_path = "{}/{}".format(dirpath, genome_file)
            model_name ="{}.xml".format(genome_file.split(".")[0])
            print(model_name)
    
            
            print("Starting Carving process")
            
            '''
            <<<Carving process by CarveMe>>>
            D. Machado et al, "Fast automated reconstruction of genome-scale metabolic models for microbial species and communities", Nucleic Acids Research, gky537, 2018. doi: https://doi.org/10.1093/nar/gky537
            '''
            universe = argv[2]
            medium = argv[3]
            initiate = argv[4]
            print_carv_status()
            run_carveme(genome_path, model_name, universe, medium, initiate)
            
            '''
            <<<Model manipulation by Cobrapy>>>
            Ebrahim, A., Lerman, J.A., Palsson, B.O. and Hyduke, D.R., 2013. COBRApy: constraints-based reconstruction and analysis for python. BMC systems biology, 7(1), pp.1-6.

            '''

            #Give the name of the solver you want to use e.g: cplex
            solver_alg = "cplex"
            cobra_config = cobra.Configuration()

            '''
            Model creation
            '''
            model =read_model(model_name)
            
            #Get basic info
            basic_info(model)

            #Checking the bounds
            print("The default flux bounds are:{}".format(cobra_config.bounds))

            #CHECK THE CONTENT OF THE MODEL

            #Reactions: Use True to print, or false to store the output
            reactions=inspect_model_reactions(model,False)

            #Metabolites: Use True to print, or false to store the output
            metabolites = inspect_model_metabolites(model,False)

            #Genes: Use True to print, or false to store the output
            genes = inspect_model_genes(model,False)

            '''
            BOF status
            
            '''
            print(model.objective)
            print(model.objective_direction)
            print("objective expression", model.objective.expression)


            '''Run FBA'''

            result = model.optimize()
            print(result)
            print(model.summary())

            '''Export sbml '''
            write_json(model, model_name)
            export_sbml(model_name, model)

            '''
            <<<Predefined maps generation by Escher>>>
            King, Z.A., Dräger, A., Ebrahim, A., Sonnenschein, N., Lewis, N.E. and Palsson, B.O., 2015. Escher: a web application for building, sharing, and embedding data-rich visualizations of biological pathways. PLoS computational biology, 11(8), p.e1004321.
            '''   
            print_generating_metabo_maps_status()
            maps_repo = "{}/all_metabo_paths".format(dirpath)
            all_maps_path = "{}/metabolism_maps".format(dirpath)
            create_dir(all_maps_path)

            maps = list_files(maps_repo)

            json_models = [x for x in list_files(dirpath) if ".json" in x]
            json_model = str([x for x in json_models if "formatted_" not in x])
            fixed_json = fix_json(json_model)
            json_model_path = "{}/{}".format(dirpath, fixed_json)

            os.chdir(all_maps_path)

            for map in maps:
                if ".json" in map:
                    map_path = "{}/{}".format(maps_repo, map)
                    escher_build(map_path, json_model_path)
            print("<<< All escher maps are in : {}>>>".format(all_maps_path))        
            os.chdir(dirpath)
            
            '''
            <<<Generation of Cytoscape compatible maps with PSAMMr>>>
            Steffensen, J.L., Dufault-Thompson, K. and Zhang, Y., 2016. PSAMM: a portable system for the analysis of metabolic models. PLoS computational biology, 12(2), p.e1004732.
            '''  

            print_cyto_status()
            model_name_f ="{}.sbml".format(model_name.split(".")[0])
            model_dir_name = "psamm_{}".format(model_name.split(".")[0])
            import_sbml_to_psamm(model_name, model_dir_name)

            os.chdir("{}/{}".format(dirpath, model_dir_name))
            psamm_path = os.getcwd()

            ''''Build data to plot'''
            generate_cyto_data(psamm_path)
    except:
        pass


def read_model(data):
    """
    This function reads a sbml file and generates the model for cobrapy

    :param data: is a file specified in the whole path under the main function
    :return: an object that can be read by the cobra suit
    """
    out_name = data.split("/")[-1][0:-5]
    if os.path.exists(out_name):
        print("{} already exists".format(out_name))
    else:
        #model = cobra.io.sbml.read_sbml_model(data, number=float, f_replace={})
        model = cobra.io.read_sbml_model(data, f_replace={})
    return(model)


def basic_info(model):
    """
    Prinst some basic info of the model

    :param model: is the model object created by crobrapy
    :return:the number o reactions, metabolites and genes in the given model
    """
    print("Basic info of the model")
    print("-----------------------")
    print("model has {} reactions".format(len(model.reactions)))
    print("model has {} metabolites".format(len(model.metabolites)))
    print("model has {} genes".format(len(model.genes)))


def inspect_model_reactions(model,status={}):
    """
    Inspect the reactions of the model

    :param model: is the model invoked
    :param status: is an argument that you must provide to the function:
        True: If wou want to print the output
        False: If you want to store the output
    :return: a dictionary containing the id of an reaction and its reaction
    """
    reactions_dict = {}
    if status == True:
        print("Reactions")
        print("----------")
        for info in model.reactions:
            print("ID:\t{}\tREACTION:\t{}".format(info.id,info.reaction))
    elif status == False:
        for info in model.reactions:
            reac_id = info.id
            reactions_dict[reac_id] = info.reaction
    return reactions_dict


def inspect_model_metabolites(model, status={}):
    """
    It inspects metabolites of the model

    :param model: is the model invoked
    :param status: Is to define what you want from the output
        True: Will print the output
        False: Will store the output in a variable
    :return: Is a dictionary containing the ID and formula of the metabolite.
    """
    metabolites_dict={}
    if status==True:
        print("Metabolites")
        print("----------")
        for info in model.metabolites:
            print("ID:\t{}\tFORMULA:\t{}".format(info.id, info.formula))
    elif status==False:
        for info in model.metabolites:
            met_id = info.id
            metabolites_dict[met_id] = info.formula
        return(metabolites_dict)


def remove_reactions(model):
    '''
    This function select and remove reactions from the model

    :param model: is cobra model
    :return:
    '''
    for i in range(len(model.metabolites)):
        m_info = model.metabolites[i]
        if "RNA" in m_info.id:
            reac = model_reactions(m_info.id)
            for i in [i.id for i in eval(reac)]:
                reaction = model.reactions.get_by_id(i)
                model.remove_reactions([reaction])


def model_reactions(met_id):
    '''
    This function removes reactions associated with selected metabolites from tge model
    :param met_id: metabolite id
    :return:
    '''
    metabolite_ass_reactions = "model.metabolites.{}.reactions".format(met_id)
    return metabolite_ass_reactions

def model_remove_metabolites(model, met_id):
    '''
    This function remove metabolites selected by remove_metabolites from the model
    :param model: is the cobra model
    :param met_id: is the metabolite id
    :return:
    '''
    for i in met_id:
        met_to_remove = "model.metabolites.{}.remove_from_model()".format(i)
        eval(met_to_remove)
def n2w(number):
    '''
    This function reads a string of numbers and convert it to a int
    :param number: is a string
    :return: converted numbers
    '''
    num2words = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', \
                 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}
    record = []
    new_number = number.lower()
    final = new_number.split(" ")
    for n in final:
        for k, v in num2words.items():
            if n == v:
                value = k
                record.append(value)
    final_value = "".join(map(str, record))
    return(final_value)


def inspect_model_genes(model,status={}):
    """
       It inspects genes of the model and check in which reactions they play role

       :param model: is the model invoked
       :param status: Is to define what you want from the output
           True: Will print the output
           False: Will store the output in a variable
       :return: Is a dictionary containing the ID and reactions associated with the gene
       """
    genes_dict = {}
    if status==True:
        print("Genes")
        print("----------")
        for info in model.genes:
            associated_ids = (i.id for i in info.reactions)
            reactions_ids = (", ".join(map(str,associated_ids)))
            print("{0}\tis associated with reactions:\t{1}".format(info.id, "{"+reactions_ids+"}"))

    elif status==False:
        for info in model.genes:
            associated_ids = (i.id for i in info.reactions)
            reactions_ids = (", ".join(map(str,associated_ids)))
            gene_id_str = info.id
            gene_id_processing = list(filter(None,gene_id_str.split('_')))[2:]
            string_id = " ".join(map(str,gene_id_processing))
            converted_str = n2w(string_id)
            gene_tag_id = "PP_{}".format(converted_str)
            genes_dict[gene_tag_id] = reactions_ids
        return(genes_dict)


def remove_metabolites(model):
    '''

    This function was inplemented to remove some metabolites from the model if needed.
    :param model: cobra model
    :return:
    '''
    met_to_remove = []
    for i in range(len(model.metabolites)):
        m_info = model.metabolites[i]
        if "RNA" in m_info.id:
            met_to_remove.append(m_info.id)
    return met_to_remove


def list_files(path):
    '''
       This function lists all file in a the path
       :param path: path of a dir
       :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path, file))]
    return (files)

def run_carveme(genome, model, universe="gramneg", medium="M9", initiate="M9"):
    '''
    This function uses carveme to build a model that can be used for manipulation

    :param genome: is the genome file that will be used to generate the model
    :param model: is the model name
    :param universe: specifies if is gram + or - organism
    :return:
    '''
    if os.path.exists(model):
        print("{} already exists".format(model))
    else:
        cmd_carve = "carve {} -u {} " \
                    "--gapfill {} -i {} --fbc2 -o {}".format(genome,universe, medium, initiate, model)
        exit_message = subprocess.check_call(cmd_carve, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} SBML model was created".format(model))

def export_sbml(model_name, model):
    '''
    This function export a cobra model to sbml file

    :param model_name: is the name that the sbml file will use
    :param model: is the cobra model
    :return:
    '''
    out_name = "{}.sbml".format(model_name.split(".")[0])
    if os.path.exists(out_name):
        print("Sbml was already generated, check:{}".format(out_name))
    else:
        cobra.io.sbml.write_sbml_model(model, out_name, f_replace={})
        print("{} was generated".format(out_name))

def write_json(model, out_name):
    '''

    This function converts cobra to json model
    :param model: is the cobra model
    :param out_name: is the json model name
    :return: the model in json format
    '''
    new_outname ="{}.json".format(out_name.split(".")[0])
    if os.path.exists(new_outname):
        print("{} already exists, moving to the next step".format(new_outname))
    else:
        cobra.io.json.save_json_model(model, new_outname, sort=False, pretty= False)
        print("{} was created.".format(new_outname))

def import_sbml_to_psamm(model,out_name):
    '''

    Imports a sbml file to psamm
    :param model: is the sbml file that you want to import
    :param out_name: the name of a dir where the info will be located
    :return: a dir with .yaml file
    '''
    if os.path.exists(out_name):
        print("{} dir already exists".format(out_name))
    else:
        cmd_psamm ="psamm-import sbml --source {} " \
                   "--dest {}/".format(model, out_name)
        exit_message = subprocess.check_call(cmd_psamm, shell = True)
        print("Exit status: {0}".format(exit_message))
        print("{0} was created".format(out_name))

def generate_cyto_data(psamm_path):
    '''
    This function creates all files needed to build a network image
    :param psamm_path: path to psamm model dir
    :return:
    '''
    files = list_files(psamm_path)
    if os.path.exists("reactions.dot"):
        print("The graph data was already generated, check tsv files for cytoscape ")
    else:
        if argv[2] != "grampos":
            cmd_vis = "psamm-model vis"
            exit_message = subprocess.check_call(cmd_vis, shell=True)
            print("The status is {}".format(exit_message))
            print("The .tsv files for cytoscape were generated")
        elif argv[2] == "grampos":
            cmd_vis_pos = "psamm-model vis --method no-fpp"
            exit_message = subprocess.check_call(cmd_vis_pos, shell=True)
            print("The status is {}".format(exit_message))
            print("The .tsv files for cytoscape were generated")


def fix_json(json_file):
    '''
    This function correct the ids of reactions and metabolites in .json model

    :param json_file: is the .json model created by cobrapy
    :return: a fixed .json file
    '''
    json_file = json_file[2:-2]
    out_file = "formatted_{}".format(json_file)
    if os.path.exists(out_file):
        print("{} already exists".format(out_file))
    else:
        with open(json_file) as file:
            text = file.read()
            text = text.replace("R_", "")
            text = text.replace("M_", "")

        original = sys.stdout
        sys.stdout = open(out_file, 'w')
        print(text)
        sys.stdout = original
    return out_file


def create_dir(path):
    '''
    This function creates new dirs
    :param path: path to the new dir
    :return:
    '''

    if os.path.exists(path):
        print("Dir {} already created".format(path))
    else:
        os.mkdir(path)


def escher_build(map_path, json_model_path):
    '''
    This function generates plots of metabolic networks

    This function uses .json model and pre-build maps to generate plots of the
    metabolic network, all files are built into the same dir.

    :param map_path: path to the map.jon
    :param json_model_path: path to .json model
    :return:
    '''
    out_path = map_path.split("/")[-1].split(".")[0]
    out_id = json_model_path.split("/")[-1].split(".")[0]
    out_name = "{}_{}.html".format(out_id, out_path)
    if os.path.exists(out_name):
        print("{} already existis".format(out_name))
    else:

        escher.rc['never_ask_before_quit'] = True

        builder = Builder(
            map_json=map_path,
            model_json=json_model_path,
        )
        builder.highlight_missing = True
        builder.save_html(out_name)
        print("{} was created".format(out_name))

def print_carv_status():
    print("---------------------------------------------------")
    print("Building Model")
    print("---------------------------------------------------")

def print_generating_metabo_maps_status():
    print("---------------------------------------------------")
    print("Generating Escher Metabolic Maps")
    print("---------------------------------------------------")

def print_cyto_status():
    print("---------------------------------------------------")
    print("Creating Cytoscape Maps")
    print("---------------------------------------------------")


def main():

        try:
            input_data = argv[1]
            if input_data == "-h" or input_data =="--help":
                control(input_data)
            if ".faa" in input_data:
                control(input_data)
        except:
            print("No input, use -h for help")





# main
if __name__ == '__main__':
    main()
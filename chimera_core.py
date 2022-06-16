# -*- coding: utf-8 -*-
# @Author: gustavotamascohotmail.com
# @Date:   2021-10-04 14:05:50
# @Last Modified by:   tamascogustavo
# @Last Modified time: 2021-10-20 10:25:17


# import statements
from tkinter.tix import Tree
from tools.path_harvest import PathColector
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
import argparse
import pandas as pd
import numpy as np
import plotly.express as px


# functions and classes
def menu():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="command")

    # core module
    core_parser = subparser.add_parser("core", help="Runs the core module")
    core_parser.add_argument("--organism", required=True,
                             help="path to the *.faa file")
    core_parser.add_argument(
        "--type", required=True, help="type of organism, options are gramneg, grampos, bacteria")
    core_parser.add_argument("--media", required=True, help="Growth media ID")

    # silencing
    core_parser = subparser.add_parser(
        "silencing", help="Runs the knockout module")
    core_parser.add_argument("--i", required=True,
                             help="path to the *.sbml model file")
    core_parser.add_argument("--type", required=True,
                             help="type of knockout target, gene or reaction")
    core_parser.add_argument("--targets", required=True,
                             help="path to csv file containing targets, one by line")

    # add new media to carveme database
    core_parser = subparser.add_parser(
        "new_media", help="Adds new media to carveme media database")
    core_parser.add_argument(
        "--table", required=True, help="path to the csv file containg the media composition")

    args = parser.parse_args()

    return args


def run_core_module(args):
    initial_path = os.getcwd()
    universe = args.type
    medium = args.media
    initiate = args.media
    genome_path = args.organism
    model_name = f'{args.organism.split("/")[-1].split(".")[0]}.sbml'

    print_carv_status()
    run_carveme(genome_path, model_name, universe, medium, initiate)

    # Start cobrapy eval
    # Give the name of the solver you want to use e.g: cplex
    solver_alg = "cplex"
    cobra_config = cobra.Configuration()

    model = read_model(model_name)

    # Get basic info
    basic_info(model)

    # Checking the bounds
    print("The default flux bounds are:{}".format(cobra_config.bounds))

    # Perform BOF checking

    print(model.objective)
    print(model.objective_direction)
    print("objective expression", model.objective.expression)

    '''Run FBA'''

    result = model.optimize()
    print(result)
    table_info = model.summary()
    print(model.summary())
    table_info.to_frame().to_csv(
        f"{model_name.split('.')[0]}_uptake_secretion_rates.csv")

    '''Export to json '''
    write_json(model, model_name)

    '''Start Escher maps production'''
    print_generating_metabo_maps_status()
    script_path = __file__
    predefined_maps = f'{"/".join(script_path.split("/")[0:-1])}/dependencies/all_metabo_paths'
    all_maps_path = "{}/metabolism_maps".format(initial_path)
    create_dir(all_maps_path)

    '''Fix json files to match Escher formatting'''
    maps = list_files(predefined_maps)

    json_models = [x for x in list_files(initial_path) if ".json" in x]
    json_model = str([x for x in json_models if "formatted_" not in x])
    fixed_json = fix_json(json_model)
    json_model_path = "{}/{}".format(initial_path, fixed_json)

    os.chdir(all_maps_path)

    for map in maps:
        if ".json" in map:
            map_path = "{}/{}".format(predefined_maps, map)
            escher_build(map_path, json_model_path)
    print("<<< All escher maps are in : {}>>>".format(all_maps_path))
    os.chdir(initial_path)

    '''Initiating psamm analysis'''

    print_cyto_status()
    model_dir_name = "psamm_{}".format(model_name.split(".")[0])
    print(model_dir_name)
    import_sbml_to_psamm(model_name, model_dir_name)
    os.chdir(model_dir_name)
    psamm_path = os.getcwd()

    ''''Build data to plot'''
    generate_cyto_data(psamm_path)

    os.chdir(initial_path)

    ''' Pathway harvest'''
    path_csv_file = f"{model_name.split('.')[0]}_enriched_paths.csv"
    if os.path.exists(path_csv_file):
        print(f"Metabolic data already collected, check: {path_csv_file}")
    else:
        model_colector = PathColector(f"{initial_path}/{model_name}")
        model_colector.parse_model_info()
        model_colector.enrich_intel()
        model_metaboism = flatten(model_colector.pathways)
        create_df(model_metaboism, f"{model_name.split('.')[0]}")
    make_plot(path_csv_file)

    exit()


def run_carveme(genome, model, universe, medium, initiate):
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
                    "--gapfill {} -i {} --fbc2 -o {}".format(
                        genome, universe, medium, initiate, model)
        exit_message = subprocess.check_call(cmd_carve, shell=True)
        print("Exit status: {0}".format(exit_message))
        print("{0} SBML model was created".format(model))


def run_knockout_module():
    pass


def add_new_media():
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


def list_files(path):
    '''
       This function lists all file in a the path
       :param path: path of a dir
       :return: all files path in a list
    '''
    files = [file for file in listdir(path) if isfile(join(path, file))]
    return (files)


def write_json(model, out_name):
    '''

    This function converts cobra to json model
    :param model: is the cobra model
    :param out_name: is the json model name
    :return: the model in json format
    '''
    new_outname = "{}.json".format(out_name.split(".")[0])
    if os.path.exists(new_outname):
        print("{} already exists, moving to the next step".format(new_outname))
    else:
        cobra.io.json.save_json_model(
            model, new_outname, sort=False, pretty=False)
        print("{} was created.".format(new_outname))


def import_sbml_to_psamm(model, out_name):
    '''

    Imports a sbml file to psamm
    :param model: is the sbml file that you want to import
    :param out_name: the name of a dir where the info will be located
    :return: a dir with .yaml file
    '''
    if os.path.exists(out_name):
        print("{} dir already exists".format(out_name))
    else:
        cmd_psamm = "psamm-import sbml --source {} " \
            "--dest {}/".format(model, out_name)
        exit_message = subprocess.check_call(cmd_psamm, shell=True)
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


def flatten(xss):
    """Makes a list of list into a single list

    Args:
        xss (list): list of lists

    Returns:
        list: return a flat list
    """
    return [x for xs in xss for x in xs]


def create_df(model_paths, name):
    """This function normalize the content of the dataframe of

    Args:
        model_paths (str): path to compound files
        name (str): name of the organism

    Returns:
        dataframe: dataframe col1 = path col2 = count of ocurrence
    """
    outname = f"{name}_enriched_paths.csv"
    model_paths = Counter(model_paths)
    model_paths = pd.DataFrame(
        model_paths.items(), columns=["Pathway", "Counts"])
    model_paths.to_csv(outname, index=False)
    return model_paths


def make_plot(df):
    """This function use a concat df to plot a clustermap

    Args:
        df (dataframe): dataframe containing the pathways and the organisms
    """
    html_outname = f"{df.split('.')[0]}_top30_pathwats.html"
    if os.path.exists(html_outname):
        pass
    else:
        combined_models_df = pd.read_csv(df, index_col=0)
        combined_models_df.index = [
            " ".join(x.split()[1:]) for x in combined_models_df.index]
        combined_models_df.drop("Metabolic pathways", axis=0, inplace=True)
        combined_models_df = combined_models_df.sort_values(
            "Counts", ascending=False)
        df = combined_models_df
        fig = px.bar(df[0:30], x='Counts', y=df[0:30].index, color=df[0:30].index, labels=dict(
            Counts="Number of compounds associated to the pathway", y="Metabolic Pathways", color="Legend"))
        fig.update_layout(
            title="Top 30 detected pathways associated to model compounds composition", title_x=0.5)

        fig.write_html(html_outname)

    exit()


def main():
    '''
    Setting the variable
    '''

    my_path = os.getcwd()
    escher.rc['never_ask_before_quit'] = True

    # Start process
    options = menu()

    CONTROLER = {"core": run_core_module,
                 "silencing": run_knockout_module,
                 "new_media": add_new_media
                 }

    execute_module = CONTROLER.get(options.command)
    execute_module(options)
    exit()

    """Special thanks to:

    <<<Predefined maps generation by Escher>>>
    King, Z.A., Dräger, A., Ebrahim, A., Sonnenschein, N., Lewis, N.E. and Palsson, B.O., 2015. Escher: a web application for building, sharing, and embedding data-rich visualizations of biological pathways. PLoS computational biology, 11(8), p.e1004321.
    
    <<<Carving process by CarveMe>>>
    D. Machado et al, "Fast automated reconstruction of genome-scale metabolic models for microbial species and communities", Nucleic Acids Research, gky537, 2018. doi: https://doi.org/10.1093/nar/gky537
   
    <<<Model manipulation by Cobrapy>>>
    Ebrahim, A., Lerman, J.A., Palsson, B.O. and Hyduke, D.R., 2013. COBRApy: constraints-based reconstruction and analysis for python. BMC systems biology, 7(1), pp.1-6.

    <<<Generation of Cytoscape compatible maps with PSAMMr>>>
    Steffensen, J.L., Dufault-Thompson, K. and Zhang, Y., 2016. PSAMM: a portable system for the analysis of metabolic models. PLoS computational biology, 12(2), p.e1004732.
    """


# main
if __name__ == '__main__':
    main()

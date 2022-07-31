#!/usr/bin/env python

'''
# @Author: gustavotamascohotmail.com
# @Date:   2021-10-04 14:05:50
# @Last Modified by:   tamascogustavo
# @Last Modified time: 2021-10-20 10:25:17

'''

# import statements
from contextlib import redirect_stdout
import json
import re
from bs4 import BeautifulSoup
from tools.path_harvest import PathColector
from tools.silencer import Silencer
import os.path
import subprocess
import sys
from itertools import islice
import cobra
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import gapfill
import os
import cobra.manipulation
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
from itertools import combinations
import ray
import requests
requests.adapters.DEFAULT_RETRIES = 5


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
    core_parser.add_argument(
        "--mediadb", help="tsv file with a new media composition for new experimental conditions. If provided you must pass the media id in --media parameter.")

    # silencing
    core_parser = subparser.add_parser(
        "silencing", help="Runs the knockout module")
    core_parser.add_argument("--i", required=True,
                             help="path to the *.sbml model file")
    core_parser.add_argument("--type", required=True,
                             help="type of knockout target, gene or reaction")
    core_parser.add_argument("--targets", required=True,
                             help="path to csv file containing targets, one by line")
    core_parser.add_argument(
        "--faa", required=True, help="path to the faa file")
    core_parser.add_argument("--mode", required=True,
                             help="Type of knockout: single or douple or all. For  double all combinations of targets will be performed")

    # complete cytoscape maps
    core_parser = subparser.add_parser(
        "harvest_path_cytoscape", help="Adds kegg pathway information to cytoscape maps")
    core_parser.add_argument(
        "--table", required=True, help="Path to reactions.edges.tsv file inside psam_* folder")
    core_parser.add_argument("--model", required=True,
                             help="path to sbml model file")

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
    run_carveme(genome_path, model_name, universe, medium, initiate, args)

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
    maps = [x for x in maps if ".py" not in x]

    #json_models = [x for x in list_files(initial_path) if ".json" in x]
    #json_model = str([x for x in json_models if "formatted_" not in x])
    json_model = f"{model_name.split('.')[0]}.json"
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
    generate_cyto_data(psamm_path, args)

    os.chdir(initial_path)

    ''' Pathway harvest'''
    print("Harvesting metabolic pathways from model. This may take a while...")
    path_csv_file = f"{model_name.split('.')[0]}_enriched_paths.csv"
    if os.path.exists(path_csv_file):
        print(f"Metabolic data already collected, check: {path_csv_file}")
    else:
        model_colector = PathColector(f"{initial_path}/{model_name}")
        model_colector.parse_model_info()
        model_colector.enrich_intel()
        model_metaboism = flatten(filter(None, model_colector.pathways))
        create_df(model_metaboism, f"{model_name.split('.')[0]}")
    make_plot(path_csv_file)


def run_carveme(genome, model, universe, medium, initiate, args):
    '''
    This function uses carveme to build a model that can be used for manipulation

    :param genome: is the genome file that will be used to generate the model
    :param model: is the model name
    :param universe: specifies if is gram + or - organism
    :return:
    '''
    if args.mediadb:
        if os.path.exists(model):
            print("{} already exists".format(model))
        else:
            cmd_carve = f"carve {genome} -u {universe} --mediadb {args.mediadb} --gapfill {medium} --fbc2 -o {model}"
            print(cmd_carve)
            exit_message = subprocess.check_call(cmd_carve, shell=True)
            print("Exit status: {0}".format(exit_message))
            print("{0} SBML model was created".format(model))
    else:
        if os.path.exists(model):
            print("{} already exists".format(model))
        else:
            cmd_carve = "carve {} -u {} " \
                "--gapfill {} -i {} --fbc2 -o {}".format(
                        genome, universe, medium, initiate, model)
            exit_message = subprocess.check_call(cmd_carve, shell=True)
            print("Exit status: {0}".format(exit_message))
            print("{0} SBML model was created".format(model))


def run_knockout_module(args):
    """This function controls the silencing module

    Args:
        args (options): are the arguments for the functions
    """
    model = args.i
    operation_to_perform = args.type
    mode = args.mode
    targets = args.targets
    faa_info = args.faa

    operator = f"{operation_to_perform}+{mode}"

    CONTROL_OPERATOR = {
        "gene+single": single_gene_knockout,
        "reaction+single": single_reaction_knockout,
        "gene+double": double_gene_knockout,
        "reaction+double": double_reaction_knockout,
        "gene+all": gene_essentiality,
        "reaction+all": reaction_essentiality
    }

    executor = CONTROL_OPERATOR.get(operator)
    executor(model, operation_to_perform, mode, targets, faa_info)


def get_silencing_targets(file):
    """Parse a txt file with gene targets       

    Args:
        file (str ): txt file with gene name, one by line
    """
    all_targets = []
    with open(file) as f:
        for line in f:
            target = line.strip()
            all_targets.append(target)
    return all_targets


def single_gene_knockout(model, operation_to_perform, mode, targets, faa_info):
    """Perform single gene knockout

    Args:
        model (str): name of the model
        operation_to_perform (str ): if gene or reaction target 
        mode (str): single or double deletion
        targets (str): txt file with targets
        faa_info (str): path to faa file
    """
    targets_to_silence = get_silencing_targets(targets)
    model_to_silence = Silencer(model, faa_info, targets_to_silence)
    sequences = model_to_silence.parse_faa_info()

    # get model genes id
    knockout_ids = model_to_silence.mine_model_gene_id(
        targets_to_silence, sequences)
    result = model_to_silence.single_gene_silence(knockout_ids)


def single_reaction_knockout(model, operation_to_perform, mode, targets, faa_info):
    """Perform single reaction knockout

    Args:
        model (str): name of the model
        operation_to_perform (str ): if gene or reaction target 
        mode (str): single or double deletion
        targets (str): txt file with targets
        faa_info (str): path to faa file
    """
    targets_to_silence = get_silencing_targets(targets)
    model_to_silence = Silencer(model, faa_info, targets_to_silence)
    result = model_to_silence.single_reaction_silence(targets_to_silence)


def double_gene_knockout(model, operation_to_perform, mode, targets, faa_info):
    """Perform double gene knockout

    Args:
        model (str): name of the model
        operation_to_perform (str ): if gene or reaction target 
        mode (str): single or double deletion
        targets (str): txt file with targets
        faa_info (str): path to faa file
    """
    targets_to_silence = get_silencing_targets(targets)
    model_to_silence = Silencer(model, faa_info, targets_to_silence)
    sequences = model_to_silence.parse_faa_info()

    # get model genes id
    knockout_ids = model_to_silence.mine_model_gene_id(
        targets_to_silence, sequences)
    all_combinations = list(combinations(list(knockout_ids.keys()), 2))
    #all_combinations = [set(x) for x in all_combinations]
    result = model_to_silence.double_gene_silence(
        knockout_ids, all_combinations)


def double_reaction_knockout(model, operation_to_perform, mode, targets, faa_info):
    """Perform double reaction knockout

    Args:
        model (str): name of the model
        operation_to_perform (str ): if gene or reaction target 
        mode (str): single or double deletion
        targets (str): txt file with targets
        faa_info (str): path to faa file
    """
    targets_to_silence = get_silencing_targets(targets)
    model_to_silence = Silencer(model, faa_info, targets_to_silence)

    all_combinations = list(combinations(targets_to_silence, 2))

    result = model_to_silence.double_reaction_silence(all_combinations)


def gene_essentiality(model, operation_to_perform, mode, targets, faa_info):
    """Perform single gene knockout for all genes in the model

    Args:
        model (str): name of the model
        operation_to_perform (str ): if gene or reaction target 
        mode (str): single or double deletion
        targets (str): txt file with targets
        faa_info (str): path to faa file
    """
    targets_to_silence = None
    model_to_silence = Silencer(model, faa_info, targets_to_silence)
    sequences = model_to_silence.parse_faa_info()

    # get model genes id
    result = model_to_silence.single_gene_silence(sequences)


def reaction_essentiality(model, operation_to_perform, mode, targets, faa_info):
    """Perform single reaction knockout for all reactions in the model

    Args:
        model (str): name of the model
        operation_to_perform (str ): if gene or reaction target 
        mode (str): single or double deletion
        targets (str): txt file with targets
        faa_info (str): path to faa file
    """
    targets_to_silence = get_silencing_targets(targets)
    model_to_silence = Silencer(model, faa_info, targets_to_silence)
    reactions_df = model_to_silence.silence_all_reactions()
    reactions_df.to_csv("reactions_essentiality.csv")


def run_harvest(args):
    print("Process initiated, this can take some time!")
    if os.path.exists("kegg_paths.json"):
        print("We collected the kegg pathways, loading info")
        with open('kegg_paths.json', 'r') as fp:
            data = json.load(fp)
        compounds_info = parse_model_info(args)
    else:
        kegg_paths = {}
        compounds_info = parse_model_info(args)
        cpd_path = []
        #teste = dict(take(10, compounds_info.items()))
        unique_kegg_set_cpd = set(list(compounds_info.values()))

        for kegg in unique_kegg_set_cpd:
            path = get_compound_info.remote(kegg)
            cpd_path.append(path)
        path_feature = ray.get(cpd_path)

        print(len(path_feature))
        clean_path_features = clean_kegg_info(path_feature)

        with open("kegg_paths.json", "w") as json_file:
            json.dump(clean_path_features, json_file, indent=4)

        with open('kegg_paths.json', 'r') as fp:
            data = json.load(fp)

    kegg_and_path = data
    bigg_and_kegg = compounds_info
    #inverted_bigg_and_kegg = {v: k for k, v in bigg_and_kegg.items()}
    add_intel_to_psamm_maps(args, kegg_and_path, bigg_and_kegg)


def clean_kegg_info(kegg_info):
    clean_info = {}
    for item in kegg_info:
        kegg_cpd = item[0]
        paths = item[1]
        if "No detection for" not in paths:
            if type(paths) != str:
                clean_info[kegg_cpd] = " | ".join(paths)
            else:
                clean_info[kegg_cpd] = paths
    return clean_info


def parse_model_info(args):
    model = cobra.io.read_sbml_model(args.model)
    compunds_intel = {}

    for cpd in model.metabolites:
        cpd_id = cpd.id
        info = cpd.annotation
        kegg_info = info.get("kegg.compound")
        if kegg_info:
            if type(kegg_info) != str:
                compunds_intel[cpd_id] = kegg_info[0]
            else:
                compunds_intel[cpd_id] = kegg_info
    return compunds_intel


@ray.remote
def get_compound_info(cpd):
    """This fuction make request to KEGG API

    Args:
        cpd (str): kegg compound ids

    Returns:
        dict: dict of pathways associated to the compound
    """

    link = "https://rest.kegg.jp/get/"
    full_path = f"{link}{cpd}"
    try:
        page = requests.get(full_path, timeout=(15, 15))
        soup = BeautifulSoup(page.content, "html.parser")
        pathway = parse_paths_bs4(soup, cpd)
        intel = (cpd, pathway)
    except (requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
        pass
        #pathway = f"No connection with {cpd}"
        #intel = (cpd, pathway)

    return intel


def parse_paths_bs4(soup_obj, cpd):
    """This function parses the request response to

    Args:
        soup_obj (soup_obj): object that contain the information of the requesr

    Returns:
        list: list of mapped pathways
    """
    paths = str(soup_obj.get_text)
    pattern = "map.*"
    result = re.findall(pattern, paths)
    # result = re.search(pattern, paths)
    if result:
        if type(result) != str:
            filter(None, result)
    else:
        result = f"No detection for {cpd}"
    return result


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
        # model = cobra.io.sbml.read_sbml_model(data, number=float, f_replace={})
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


def generate_cyto_data(psamm_path, args):
    '''
    This function creates all files needed to build a network image
    :param psamm_path: path to psamm model dir
    :return:
    '''
    files = list_files(psamm_path)
    if os.path.exists("reactions.dot"):
        print("The graph data was already generated, check tsv files for cytoscape ")
    else:
        if args.type != "grampos":
            cmd_vis = "psamm-model vis"
            exit_message = subprocess.check_call(cmd_vis, shell=True)
            print("The status is {}".format(exit_message))
            print("The .tsv files for cytoscape were generated")
        elif args.type == "grampos":
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
    json_file = json_file
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
    """This function plots the top30 pathways detected in the model

    Args:
        df (dataframe): dataframe containing the pathways in the model
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


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def add_intel_to_psamm_maps(args, path_info, bigg_info):
    """This function will add pathway information to psamm edges files

    Args:
        args (class): holds all options privided by user.
        path_info (dict): k: kegg id , v : pathway
        bigg_info (dict): k: bigg id , v : kegg id
    """
    if os.path.exists(f'{args.model.split(".")[0]}_reactions_edges.tsv'):
        print(f'{args.model.split(".")[0]}_reactions_edges.tsv already exists')
    else:
        with open(args.table) as file:
            file_content = parse_edges_info(args, file, path_info, bigg_info)


def parse_edges_info(args, file, path_dict, cpd_dict):
    with open(f'{args.model.split(".")[0]}_reactions_edges.tsv', 'w') as f:
        with redirect_stdout(f):
            print(
                'source\ttarget\tdirection\tpenwidth\tstyle\tsource_pathway\ttarget_pathway')
            for line_index, line in enumerate(file):
                if line_index >= 1:
                    newline = line.strip().split("\t")
                    source = newline[0]
                    target = newline[1]
                    direction = newline[2]
                    penwidth = newline[3]
                    style = newline[4]
                    pathway_source = incorporate_pathway_info(
                        source, path_dict, cpd_dict)
                    pathway_target = incorporate_pathway_info(
                        target, path_dict, cpd_dict)
                    print(
                        f'{source}\t{target}\t{direction}\t{penwidth}\t{style}\t{pathway_source}\t{pathway_target}')


def incorporate_pathway_info(source, path_dict, cpd_dict):
    source = source.replace("[c]", "").replace("[e]", "").replace("[p]", "")
    kegg_name = cpd_dict.get(source, "no match")
    associated_path = path_dict.get(kegg_name, "no pathway detected")
    return associated_path


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
                 "harvest_path_cytoscape": run_harvest
                 }

    execute_module = CONTROLER.get(options.command)
    execute_module(options)

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

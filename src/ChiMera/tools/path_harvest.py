
import cobra
import requests
from bs4 import BeautifulSoup
import re
from collections import Counter
import pandas as pd
import ray
import numpy as np

class PathColector:
    def __init__(self, model):
        self.model = model
    
    def parse_model_info(self):
        model = cobra.io.read_sbml_model(self.model)
        info_list = []

        for cpd in model.metabolites:
            info = cpd.annotation
            kegg_info = info.get("kegg.compound")
            info_list.append(kegg_info)

        info_list = list(filter(None, info_list))
        flat_list = []
        for item in info_list:
            if type(item) != str:
                item = item[0]
                flat_list.append(item)
            else:
                flat_list.append(item)

        self.metabolites = set(flat_list)

        return self.metabolites

    @classmethod
    @ray.remote
    def get_compound_info(cpd):
        """This fuction make request to KEGG API

        Args:
            cpd (str): kegg compound ids

        Returns:
            list: list of pathways associated to the compound
        """
        try:
            link = "https://rest.kegg.jp/get/"
            full_path = f"{link}{cpd}"
            page = requests.get(full_path, timeout=(15,15))
            soup = BeautifulSoup(page.content, "html.parser")
            pathway = PathColector.parse_paths(soup)
            if pathway:
                return pathway
        except (requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
            pass
    
        
    
    
    def enrich_intel(self):
        """This function collect pathay information from a model file

        Args:
            metabolites (list): list of kegg compound ids

        Returns:
            list: list of pathways
        """
        features = [PathColector.get_compound_info.remote(cpd) for cpd in self.metabolites]
        pathway_list = ray.get(features)
        self.pathways = pathway_list
        return self.pathways
        


    

    def parse_paths(soup_obj):
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
        return result



    
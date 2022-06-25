import pandas as pd
import cobra
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)


class Silencer:
    """This class control the knockout module
    """

    def __init__(self, model, composition, targets):
        self.model = cobra.io.read_sbml_model(model)
        self.fasta = composition
        self.targets = targets

    def parse_faa_info(self):
        """This  function parses a fasta file 

        Returns:
            dict: key:locus name v: gene name
        """
        with open(self.fasta) as faa_file:
            faa_info = {}
            for line in faa_file:
                if line.startswith(">"):
                    if "[gene=" in line:
                        index = line.strip().split()[0].split("|")[
                            0].strip(">")
                        gene_locus = line.strip().split()[0].split("|")[1]
                        gene_locus = "{}_{}".format(
                            index, gene_locus.replace(".", "_"))
                        gene_id = line.strip().split()[1].split("=")[1][0:-1]
                        faa_info[gene_locus] = gene_id
        return faa_info

    def mine_model_gene_id(self, targets, genes_info):
        """This function gets gene locus for the targets genes

        Args:
            targets (list): list of target genes
            genes_info (dict): dict produced by parse_faa

        Returns:
            dict: key:locus name v: gene name for the target genes
        """
        knockout_dict = {}
        for target in targets:
            if target in genes_info.values():
                knock_gene = list(genes_info.keys())[
                    list(genes_info.values()).index(target)]
                knockout_dict[knock_gene] = target
        return knockout_dict

    def single_gene_silence(self, targets):
        """This function perform single gene deletion

        Args:
            targets (dict): key:locus name v: gene name for the target genes
        """
        print('complete model growth rate: ', self.model.optimize())
        single = single_gene_deletion(self.model)
        for target in targets.keys():
            with self.model:
                try:
                    teste = eval(f"self.model.genes.{target}")
                    print(f"\n\nThe knockout effect of {targets.get(target)}")
                    print(single.knockout[teste])
                except AttributeError:
                    pass

    def double_gene_silence(self, targets, targets_combination):
        """This function perform double gene deletion

        Args:
            targets (dict of targets): key:locus name v: gene name for the target genes
            targets_combination (list): combination of all genes to silence 2 by 2
        """
        print('complete model growth rate: ', self.model.optimize())

        for target in targets_combination:
            gene1 = target[0]
            gene2 = target[1]
            with self.model:
                try:
                    silence1 = eval(f"self.model.genes.{gene1}.knock_out()")
                    silence2 = eval(f"self.model.genes.{gene2}.knock_out()")
                    print(
                        f"\n\nThe knockout effect of {targets.get(gene1)} + {targets.get(gene2)}")
                    print(self.model.optimize())
                except AttributeError:
                    pass

    def single_reaction_silence(self, targets):
        """This function perform single gene deletion

        Args:
            targets (list): all reactions targets, must be upper case
        """
        print('complete model growth rate: ', self.model.optimize())
        single = single_reaction_deletion(self.model)
        for target in targets:
            with self.model:
                try:
                    teste = eval(f"self.model.reactions.{target}")
                    print(f"\n\nThe knockout effect of {target}")
                    print(single.knockout[teste])
                except AttributeError:
                    pass

    def double_reaction_silence(self, targets_combination):
        """This function perform double reac deletion

        Args:
            targets (list): all reactions targets, must be upper case
        """
        print('complete model growth rate: ', self.model.optimize())
        for target in targets_combination:
            reac1 = target[0]
            reac2 = target[1]
            with self.model:
                try:
                    silence1 = eval(
                        f"self.model.reactions.{reac1}.knock_out()")
                    silence2 = eval(
                        f"self.model.reactions.{reac2}.knock_out()")
                    print(
                        f"\n\nThe knockout effect of {reac1} + {reac2}")
                    print(self.model.optimize())
                except AttributeError:
                    pass

    def silence_all_reactions(self):
        """Perform reactions essentiality detection

        Returns:
            dataframe: dataframe of the importance of the reaction in the organism
        """
        print('complete model growth rate: ', self.model.optimize())
        result = (single_reaction_deletion(
            self.model, self.model.reactions[0:], processes=6).round(4))
        print("Check eactions_essentiality.csv for more details")
        return result

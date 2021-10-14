# Chimera: An easy-to-use pipeline for Genome-based Metabolic Network reconstruction, evaluation, and visualization.
![](https://img.shields.io/badge/<python>-<3.7>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<carveme>-<1.4.1>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<diamond>-<v2.0.9.147>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<cobrapy>-<0.22.1>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<escher>-<1.7.3>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<psamm>-<1.1.2>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<cplex>-<studio2010>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)

![](https://imgur.com/a/AQYwl4O)

Several genome-scale metabolic reconstruction (GSMR) tools have been developed in the last decades. These tools have helped to construct many metabolic models, which had contributed to a variety of fields, e.g., genetic engineering, drug discovery, prediction of phenotypes, and other model-driven discoveries. However, the use of these tools requires a high level of programming capabilities. Multiple steps need to be accounted for, before the generation of a functional model able to produce predictions. Another limitation is the lack of a visualization module, something that can contribute to the understanding of the metabolic network. Therefore, there is a scarcity of user-friendly tools that can be used in daily routine, providing insights about the metabolic network of a target organism for researcher’s groups.

Here we present a novel tool, Chimera, which combines the most efficient tools in model reconstruction, prediction, and visualization and also implements new in-house algorithms for database integration and data manipulation. 

## Our aim with Chimera

- Produce an organism-specific model based on the CarveMe algorithm 

- Manage the model and perform growth predictions with COBRApy 

- Create visualization for the metabolic network using PSAMM and Escher 

- Addition of pathway information to metabolic maps using in house algorithm 

- Perform single and double, gene and reaction, knockout in the organism 


All of that is obtained through automation and connection of the most used tools in the literature, creating a pipeline that is easy to use, helping those researchers with none or small programming capabilities.


## Detail about the scripts

You must download the folder or the git repository. In the folder there are 3 scripts, one extra folder and 2 txt files with the dependencies to install

```
chimera_core.py: perform the model creation evaluation and visualization files

translator_using_bigg.py: is used to add KEGG pathway description in the genome-scale metabolic map

simulating_knockouts.py: is used to perform gene and reactions knockout

all_metabo_paths: a folder containing pre-defined metabolic maps for a general overview of the pathways present in the target organism
```
`conda and pip requirements can be installed following the instructions`
### Installation

Before installing you need to install CPLEX solver from IMB. Click [here](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students) to access the academic license. This documentation was created with `IBM ILOG CPLEX Optimization Studio V20.10`.

After download, follow the required system installation. For Linux, make sure to export your installation path

```
export PATH=/PATH_TO_CPLEX/cplex/bin/x86-64_linux/:$PATH
```

Install the conda env

```
conda env create -f environment.yml
# activate the environment
conda activate chimera
```

Navigate to `python cplex` the installation folder:

```
cd /PATH_TO_CPLEX/cplex/python/3.7/x86-64_linux/
```

And then install

```
python setup.py install
```
_OBS_: If you have troubles installing the API, check : https://www.ibm.com/docs/en/icos/12.8.0.0?topic=cplex-setting-up-python-api

## Update

If a new library needs to be installed, don't forget to update the environment.yml file

```
conda env export | grep -v "^prefix: " > environment.yml
```

### Usage

Before using, activate the conda env

```
conda activate chimera
```

# To access the help page:

```
python chimera_core.py -h
```

# To run the test on model  __Escherichia coli__:

```
python chimera_core.py input_examples/faa_file/e_coli_test.faa gramneg LB LB
```
# To perfom pathway annotation using KEGG as reference to the Cytoscape maps we can use:

```
python3 translator_using_bigg.py e_coli
```
This command will annotate the pathway for the compounds in the metabolic map. During the process a few warning messages can be displayed due to multiple API request. However, thats not an error.

# To perform gene or reaction knockout we can use on the model __Escherichia coli__:

**For specific genes**

Single gene deletion for all genes in the knockout_genes_list.txt:

```
python3 simulating_knockouts.py -sg input_examples/faa_file/e_coli_test.faa input_examples/reations_gene_to_delete/knockout_genes_list.txt

```
Double gene deletion for all genes in the knockout_genes_list.txt:

```
python3 simulating_knockouts.py -dg input_examples/faa_file/e_coli_test.faa input_examples/reations_gene_to_delete/knockout_genes_list.txt

```

**For specific reactions**

Single reaction deletion for all reactions in the knockout_reactions_list.txt:

```
python3 simulating_knockouts.py -sr input_examples/faa_file/e_coli_test.faa input_examples/reations_gene_to_delete/knockout_reactions_list.txt

```
Double reaction deletion for all reactions in the knockout_reactions_list.txt:

```
python3 simulating_knockouts.py -dr input_examples/faa_file/e_coli_test.faa input_examples/reations_gene_to_delete/knockout_reactions_list.txt

```
Specifc Double reaction deletion. Here we only evaluate the first 2 reactions in  knockout_reactions_list.txt:

```
python3 simulating_knockouts.py -tdr input_examples/faa_file/e_coli_test.faa input_examples/reations_gene_to_delete/knockout_reactions_list.txt

```
**Gene and reaction essenciality **

Here we evaluate the individual impact in growth due to a single gene or reaction in the model.

To evaluate gene essenciality:

```
python3 simulating_knockouts.py -cg

```
To evaluate reaction essenciality:

```
python3 simulating_knockouts.py -cr

```

**GIVE EXAMPLES OF OTHER SCRIPTS AND DISCUSS RESULTS**


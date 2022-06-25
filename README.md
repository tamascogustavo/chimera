# ChiMera: An easy-to-use pipeline for Genome-based Metabolic Network reconstruction, evaluation, and visualization.
![](https://img.shields.io/badge/<python>-<3.7>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<carveme>-<1.4.1>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<diamond>-<v2.0.9.147>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<cobrapy>-<0.22.1>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<escher>-<1.7.3>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<psamm>-<1.1.2>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<cplex>-<studio20.10>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)

![](https://imgur.com/a/AQYwl4O)

Several genome-scale metabolic reconstruction (GSMR) tools have been developed in the last decades. These tools have helped to construct many metabolic models, which had contributed to a variety of fields, e.g., genetic engineering, drug discovery, prediction of phenotypes, and other model-driven discoveries. However, the use of these tools requires a high level of programming capabilities. Multiple steps need to be accounted for, before the generation of a functional model able to produce predictions. Another limitation is the lack of a visualization module, something that can contribute to the understanding of the metabolic network. Therefore, there is a scarcity of user-friendly tools that can be used in daily routine, providing insights about the metabolic network of a target organism for researcher’s groups.

Here we present a novel tool, Chimera, which combines the most efficient tools in model reconstruction, prediction, and visualization and also implements new in-house algorithms for database integration and data manipulation. 

## Our aim with Chimera

- Produce an organism-specific model based on the CarveMe algorithm 

- Manage the model and perform growth predictions with COBRApy 

- Create visualization for the metabolic network using PSAMM and Escher 

- Add pathway information to metabolic maps using in house algorithm 

- Perform single and double, gene and reaction, knockout in the organism 


All of that is obtained through automation and connection of the most used tools in the literature, creating a pipeline that is easy to use, helping those researchers with none or small programming capabilities.



### Installation

Verify if you have anaconda installed in your machine. The tool can be downloaded [here](https://www.anaconda.com/products/individual).

Before installing you need to install CPLEX solver from IMB. Click [here](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students) to access the academic license. This documentation was created with `IBM ILOG CPLEX Optimization Studio V20.10`.

After download, follow the required system installation. For Linux, make sure to export your installation path

```
export PATH=/PATH_TO_CPLEX/cplex/bin/x86-64_linux/:$PATH
```
Obs: The command above will export the installation path of cplex to $PATH shell variable. Pay attention during the installation procedure, it will ask you where to install cplex (<PATH_TO_CPLEX>)

**Install the conda env**

There are 2 envs, use the one that suits you better:

environment.yml for mac
environment_linux.yml for linux


```
conda env create -f <environment.yml>
# activate the environment
conda activate chimera
```
Alternative intall
```
__Creating env__
conda create --name chimera python=3.7
conda activate chimera
__Install dependencies__
conda install pip
conda install -c bioconda diamond
pip install carveme
pip install cobra
pip install escher
pip install git+https://github.com/zhanglab/psamm-import.git
conda install -c conda-forge biopython
conda install -c conda-forge bioservices

```

Navigate to `python cplex` the installation folder:

```
cd /PATH_TO_CPLEX/cplex/python/3.7/x86-64_linux/
```

And then install

```
python setup.py install
```
_OBS_: If you have troubles installing the API, check: https://www.ibm.com/docs/en/icos/12.8.0.0?topic=cplex-setting-up-python-api .

_OBS_: Some users reported an error during the installation of the API.  "Unknown distribution option: zip_fafe. Could not create build: permission denied".

To solve, instead the command above, use:

```
python setup.py build -b /home_directory/
```
Or, you can follow [this tutorial](https://askubuntu.com/questions/263450/folders-showing-lock-icon).


## Update

If a new library needs to be installed, don't forget to update the environment.yml file

```
conda env export | grep -v "^prefix: " > environment.yml
```

## Usage

Before using, activate the conda env

```
conda activate chimera
```
__Intended usage__: You should have a folder with the 3 main scripts, the all_metado_paths folder and a faa file for your organism. 

### To add new media for model creation:

You can check the input_examples folder to see how to build your new_media.tsv. The ids must be BiGG Ids.

__If fail happens__:

Failed to gapfill model for medium <your_added_media> means there is no possible solution to make the organism grow on that medium.

Your may have a medium that lacks elements (iron, magnesium, zinc, etc), which are necessary for biomass formation. You can use the M9 media composition
as a template to ensure that all necessary elements are present.

You can add all the trace elements and repeat the reconstruction to check if your organism can or not grow in the provided media.

### To access the help page:

```
python chimera_core.py -h
```

### To run the test on model  __Escherichia coli__:

```
python chimera_core.py core --organism input_examples/faa_file/e_coli_test.faa --type gramneg --media M9
```
You can also use a pre-buit model, overstepping the model creation, directly producing FBA predictions and Visualization. You just need to have .xml/.sbml model in your folder, with the same prefix as your faa file. 

The same command is used in this case.
```
python chimera_core.py core --organism input_examples/faa_file/e_coli_test.faa --type gramneg --media M9
```
Due to annotation discrepancies, cytoscape compatible file may fail to be created. This is due to id mismatch of the provided model.
If you still want to perform the graph creation, inside psamm folder type the following code:

```
psamm-model vis --method no-fpp
```

### To perform pathway annotation using KEGG as a reference to the Cytoscape maps we can use:

OBS: We need the previous step to be executed before running this module.

```
chimera_core.py harvest_path_cytoscape --table <Path to reactions.edges.tsv file inside psam_* folder> --model <path to sbml model file>
```

### To perform gene or reaction knockout we can use on the model __Escherichia coli__:

OBS: Using *.faa file from Prodigal may cause inconsistencies due to the annotation labeling.

Best results are produced with annotation made with NCBI-PGAP.

```
chimera_core.py silencing --i I --type TYPE --targets TARGETS
                                 --faa FAA --mode MODE

optional arguments:
  -h, --help         show this help message and exit
  --i I              path to the *.sbml model file
  --type TYPE        type of knockout target, gene or reaction
  --targets TARGETS  path to csv file containing targets, one by line
  --faa FAA          path to the faa file
  --mode MODE        Type of knockout: single or douple or all. For double all
                     combinations of targets will be performed

```

## Build your custom maps Escher

If you want to build your own custom map based on your metabolic evidence you can do that using your formatted JSON model at [Escher website](https://escher.github.io/#/)

For instructions on how to do that, check the [tutorial](https://youtu.be/XQRbSkvMpN4).

## Build your custom maps Cytoscape

If you want to build your own custom map based on your metabolic evidence you can do that using your <org_name>reactions_edges.tsv file in the main directory of the results

For instructions on how to do that, check the [tutorial](https://youtu.be/M7SNCnPwqF0).


## Outputs 

+ In the main folder:

    + _enriched_paths.csv = csv file containing the metabolic paths of your organism
    + _enriched_paths_top30_pathways.html = plot of the top 30 detected paths
    + _uptake_secretion_rates.csv = secreteed (- values) and uptaked (+ values) compounds and their fluxes.
    + models in format json, sbml and a formated json (compatible with Escher)

+ Inside psamm folder:

    + yaml file containing compounds info
    + yaml file containing reactions info
    + yaml file containing exchange reactions info
    + yaml file containing model overall info
    + reactions.edges.tsv and reactions.nodes.tsv = files that can be loaded into Cytoscape

+ The harvest_path module will produce a *organism*_reactions_edges.tsv in the main folder.
    + This file can be loaded in Cytoscape, however now you can search for specific pathways

+ The silencing module will print to terminal the result of the knockout, only when performing knockout of all reactions in the model, a reactions_essentiality.csv will be generated

## LICENSE

The Chimera source is released under both the GPL and LGPL licenses version 3 or later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License or the GNU Lesser General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Check the [Chimera License](LICENSE). Gustavo Tamasco.

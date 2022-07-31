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


## Full tutorial

Check our `ChiMera_class` notebook !

### Installation
To install Chimera you can just check [CHiMera pip instruction](https://pypi.org/project/ChiMera-ModelBuilder/)

```
pip install ChiMera-ModelBuilder
```

_OBS_: If you plan to use envs, create one with python 3.7 or 3.6. And make sure that you install cplex in the same environment.

You also need to install CPLEX solver from IMB. Click [here](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students) to access the academic license. This documentation was created with `IBM ILOG CPLEX Optimization Studio V20.10`.

After download, follow the required system installation. For Linux, make sure to export your installation path

```
export PATH=/PATH_TO_CPLEX/cplex/bin/x86-64_linux/:$PATH
```
Obs: The command above will export the installation path of cplex to $PATH shell variable. Pay attention during the installation procedure, it will ask you where to install cplex (<PATH_TO_CPLEX>)

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


### To add new media for model creation:

You can check the input_examples folder to see how to build your new_media.tsv. The ids must be BiGG Ids.

```
python chimera_core.py core --organism input_examples/faa_file/e_coli_test.faa --type gramneg --mediadb <new_media.tsv> --media <new_media_id>
```

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

## All comand options

# Core
```
usage: chimera_core.py core [-h] --organism ORGANISM --type TYPE --media MEDIA
                            [--mediadb MEDIADB]

optional arguments:
  -h, --help           show this help message and exit
  --organism ORGANISM  path to the *.faa file
  --type TYPE          type of organism, options are gramneg, grampos,
                       bacteria
  --media MEDIA        Growth media ID
  --mediadb MEDIADB    tsv file with a new media composition for new
                       experimental conditions. If provided you must pass the
                       media id in --media parameter.
```
# Silencing
```
usage: chimera_core.py silencing [-h] --i I --type TYPE --targets TARGETS
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
# Harvest pathway

```
usage: chimera_core.py harvest_path_cytoscape [-h] --table TABLE --model MODEL

optional arguments:
  -h, --help     show this help message and exit
  --table TABLE  Path to reactions.edges.tsv file inside psam_* folder
  --model MODEL  path to sbml model file
```
## BIOMASS DEFAULT FUNCTIONS
### Gram-negative
```
Growth: 0.000223 10fthf_c + 0.513689 ala__L_c + 0.000223 amet_c + 0.295792 arg__L_c + 0.241055 asn__L_c + 0.241055 asp__L_c + 54.124831 atp_c + 0.005205 ca2_c + 0.005205 cl_c + 0.000576 coa_c + 0.0001 cobalt2_c + 0.133508 ctp_c + 0.000709 cu2_c + 0.09158 cys__L_c + 0.026166 datp_c + 0.027017 dctp_c + 0.027017 dgtp_c + 0.026166 dttp_c + 0.000223 fad_c + 0.006715 fe2_c + 0.007808 fe3_c + 0.26316 gln__L_c + 0.26316 glu__L_c + 0.612638 gly_c + 0.215096 gtp_c + 48.601527 h2o_c + 0.094738 his__L_c + 0.290529 ile__L_c + 0.195193 k_c + 0.019456 kdo2lipid4_p + 0.450531 leu__L_c + 0.343161 lys__L_c + 0.153686 met__L_c + 0.008675 mg2_c + 0.000223 mlthf_c + 0.000691 mn2_c + 0.0001 mql8_c + 0.013894 murein5px4p_p + 0.001831 nad_c + 0.000447 nadp_c + 0.017868 pe160_c + 0.045946 pe160_p + 0.054154 pe161_c + 0.02106 pe161_p + 0.185265 phe__L_c + 0.000223 pheme_c + 0.221055 pro__L_c + 0.000223 pydx5p_c + 0.000223 ribflv_c + 0.215792 ser__L_c + 0.000223 sheme_c + 0.004338 so4_c + 0.000223 thf_c + 0.000223 thmpp_c + 0.253687 thr__L_c + 0.056843 trp__L_c + 0.137896 tyr__L_c + 0.144104 utp_c + 0.423162 val__L_c + 0.000341 zn2_c --> 53.95 adp_c + 53.95 h_c + 53.945662 pi_c + 0.773903 ppi_c
```
### Gram-positive
```
Growth: 0.000223 10fthf_c + 0.513689 ala__L_c + 0.000223 amet_c + 0.295792 arg__L_c + 0.241055 asn__L_c + 0.241055 asp__L_c + 54.124831 atp_c + 0.005205 ca2_c + 0.005205 cl_c + 0.000576 coa_c + 0.0001 cobalt2_c + 0.133508 ctp_c + 0.000709 cu2_c + 0.09158 cys__L_c + 0.026166 datp_c + 0.027017 dctp_c + 0.027017 dgtp_c + 0.026166 dttp_c + 0.000223 fad_c + 0.006715 fe2_c + 0.007808 fe3_c + 0.26316 gln__L_c + 0.26316 glu__L_c + 0.612638 gly_c + 0.001 gtca1_45_BS_c + 0.001 gtca2_45_BS_c + 0.001 gtca3_45_BS_c + 0.215096 gtp_c + 48.601527 h2o_c + 0.094738 his__L_c + 0.290529 ile__L_c + 0.195193 k_c + 0.450531 leu__L_c + 5e-05 lipo1_24_BS_c + 5e-05 lipo2_24_BS_c + 5e-05 lipo3_24_BS_c + 5e-05 lipo4_24_BS_c + 0.343161 lys__L_c + 0.153686 met__L_c + 0.008675 mg2_c + 0.000223 mlthf_c + 0.000691 mn2_c + 0.0001 mql8_c + 0.001831 nad_c + 0.000447 nadp_c + 0.01 peptido_BS_c + 0.185265 phe__L_c + 0.221055 pro__L_c + 0.000223 pydx5p_c + 0.000223 ribflv_c + 0.215792 ser__L_c + 0.004338 so4_c + 0.000223 thf_c + 0.000223 thmpp_c + 0.253687 thr__L_c + 0.056843 trp__L_c + 0.137896 tyr__L_c + 0.144104 utp_c + 0.423162 val__L_c + 0.000341 zn2_c --> 53.95 adp_c + 53.95 h_c + 53.945662 pi_c + 0.773903 ppi_c
```
## LICENSE

The Chimera source is released under both the GPL and LGPL licenses version 3 or later. You may choose which license you choose to use the software under.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License or the GNU Lesser General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Check the [Chimera License](LICENSE). Gustavo Tamasco.

## Special thanks and Please cite them !


 + Predefined maps generation by Escher
 
    King, Z.A., Dräger, A., Ebrahim, A., Sonnenschein, N., Lewis, N.E. and Palsson, B.O., 2015. Escher: a web application for building, sharing, and embedding data-rich visualizations of biological pathways. PLoS computational biology, 11(8), p.e1004321.
    
+ Carving process by CarveMe

    D. Machado et al, "Fast automated reconstruction of genome-scale metabolic models for microbial species and communities", Nucleic Acids Research, gky537, 2018. doi: https://doi.org/10.1093/nar/gky537   
+ Model manipulation by Cobrapy

    Ebrahim, A., Lerman, J.A., Palsson, B.O. and Hyduke, D.R., 2013. COBRApy: constraints-based reconstruction and analysis for python. BMC systems biology, 7(1), pp.1-6.
    
+ Generation of Cytoscape compatible maps with PSAMM

    Steffensen, J.L., Dufault-Thompson, K. and Zhang, Y., 2016. PSAMM: a portable system for the analysis of metabolic models. PLoS computational bio
    

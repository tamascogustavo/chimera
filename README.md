### Chimera: An easy-to-use pipeline for Genome-based Metabolic Network reconstruction, evaluation, and visualization.
![](https://img.shields.io/badge/<python>-<3.7>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<carveme>-<1.4.1>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<diamond>-<v2.0.9.147>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<cobrapy>-<0.22.1>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<escher>-<1.7.3>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<psamm>-<1.1.2>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a) ![](https://img.shields.io/badge/<cplex>-<studio2010>-informational?style=flat&logo=<LOGO_NAME>&logoColor=white&color=2bbc8a)



Several genome-scale metabolic reconstruction (GSMR) tools have been developed in the last decades. These tools have helped to construct many metabolic models, which had contributed to a variety of fields, e.g., genetic engineering, drug discovery, prediction of phenotypes, and other model-driven discoveries. However, the use of these tools requires a high level of programming capabilities. Multiple steps need to be accounted for, before the generation of a functional model able to produce predictions. Another limitation is the lack of a visualization module, something that can contribute to the understanding of the metabolic network. Therefore, there is a scarcity of user-friendly tools that can be used in daily routine, providing insights about the metabolic network of a target organism for researcher’s groups.

Here we present a novel tool, Chimera, which combines the most efficient tools in model reconstruction, prediction, and visualization and also implements new in-house algorithms for database integration and data manipulation. 

### Our aim with Chimera

Produce an organism-specific model based on the CarveMe algorithm 

Manage the model and perform growth predictions with COBRApy 

Create visualization for the metabolic network using PSAMM and Escher 

Addition of pathway information to metabolic maps using in house algorithm 

Perform single and double, gene and reaction, knockout in the organism 


All of that is obtained through automation and connection of the most used tools in the literature, creating a pipeline that is easy to use, helping those researchers with none or small programming capabilities.


###Detail about the scripts

You must download the folder or the git repository. In the folder there are 3 scripts and one extra folder:

chimera_core.py: perform the model creation evaluation and visualization files

translator_using_bigg.py: is used to add KEGG pathway description in the genome-scale metabolic map

simulating_knockouts.py: is used to perform gene and reactions knockout

all_metabo_paths: a folder containing pre-defined metabolic maps for a general overview of the pathways present in the target organism



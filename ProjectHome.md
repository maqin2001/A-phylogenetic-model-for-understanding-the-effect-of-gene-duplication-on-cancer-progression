# **Cancer Evaluation: A Phylogenetic Model for Understanding The Effect of Gene Duplication on Cancer Progression** #

## Introduction ##
In the Bayesian phylogenetic model, the progressive stages of cancer are described as the tips of a phylogenetic tree. The duplication/deletion events may occur randomly along the branches of the tree. This model, together with duplication data, can estimate the duplication and deletion rates as cancer advances.

## How to Use This Program ##
To compile the source code, open "makefile" and change "ARCHITECTURE ?= mac" to the platform (windows, mac, unix) of your computer. Then type make.

To run the program, type './cancer' command. It will generate an output file .out in which the loglikelihood and model parameter 'm' are sampled from the MCMC algorithm. The posterior distribution of the phylogenetic tree is saved in the output file .tre.


## Contact ##
Any questions, problems, bugs are welcome and should be dumped to Qin Ma <**_qin.ma@sdstate.edu_**>

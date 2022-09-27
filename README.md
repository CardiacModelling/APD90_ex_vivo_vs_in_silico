# Benchmark of undiseased human ventricular cell models

Hi fellow scientist!

Welcome to this repo supporting the paper "Comparison of in silico predictions of action potential duration in response to inhibition of IKr and ICaL with new human experimental data", (Barral et al. 2022). This is the place where you will find all the necessary resources to reproduce our results. To make sure that everything runs smoothly, please beforehand:
  - Install the Python Package Myokit (http://myokit.org/) (Michael Clerx, Pieter Collins, Enno de Lange, and Paul GA Volders.  Myokit:  a simple interface to cardiac cellularelectrophysiology.Progress in biophysics and molecular biology, 120(1-3):100–114, 2016)
  - Install the Python Package SABS_PKPD (https://github.com/rcw5890/SABS_project)
  
Once you have the required packages, you should be able to run the scripts included in the repo, and reproduce our figures. Please be aware that the script to generate the 2D maps uses the multiprocessing library of Python, **which may not be compatible with some Windows machines.** You may need to disable the parallelisation of the evaluation of the map then. The generation of 2D maps (those for ORd based models in particular) may take a few hours so be patient ;)

You might also want to quickly run the models using OpenCOR, a simulation interface for CellML models (https://opencor.ws/) (Garny, Alan, and Peter J. Hunter. "OpenCOR: a modular and interoperable approach to computational biology." Frontiers in physiology 6 (2015): 26.)

Have a great day of modelling and simulation!

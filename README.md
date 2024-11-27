# Benchmark of undiseased human ventricular cell models

Hi fellow scientist!

Welcome to this repo supporting the paper "Comparison of in silico predictions of action potential duration in response to inhibition of IKr and ICaL with new human experimental data", (Barral et al. 2024). This is the place where you will find all the necessary resources to reproduce our results. 

## Creating the necessary Python environment
If you don't have it yet on your machine, start by installing Python. 

### Option 1

Create a python3 virtual environment, start it, and install a few dependencies:
```
python3 -m venv .venv
source .venv/bin/activate
pip install myokit pandas scipy
```
If anything fails to build/run the original dependency versions can be seen in the environment.yml mentioned in Option 2.

### Option 2
Use the Anaconda suite (https://www.anaconda.com/products/distribution).

To make sure that all scripts runs smoothly, please create first the virtual environment with necessary libraries, using the ```environment.yml``` file. To do so, browse to this folder's location and run the command:
```conda env create -f environment.yml```.

In case you have troubles with Myokit, please make sure that you satisfied the requirements (http://myokit.org/install).
  
## Run the scripts to reproduce our results
Once you have the required packages, you should be able to run the scripts included in the repo, and reproduce our figures. 
To pass the corresponding arguments to the Python scripts, please refer to the first section of each script. 

The simulation of DAPD 2-D maps was computationally expensive.
We ran our scripts on an HPC (running with LSF) to parallelise the runs, so we could submit arrays of executions of the scripts, which tremendously sped up the execution.
Note that the generation of the 2-D maps on a local machine may take a few hours  (those for ORd-similar models in particular) so be patient ;)

## Visualisation of model outputs with OpenCOR
Additionally, the models can be run easily with OpenCOR, a simulation interface for CellML models (https://opencor.ws/) (Garny, Alan, and Peter J. Hunter. "OpenCOR: a modular and interoperable approach to computational biology." Frontiers in physiology 6 (2015): 26.). 
This enables quick visualisation of model outputs.

Have a great day of modelling and simulation!

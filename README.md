# Multi-strain phage induced clearance of bacterial infections


This repository contains the source code associated with the manuscript

Marchi, Ngoc Minh, Debarbieux, Weitz: Multi-strain phage induced clearance of bacterial infections 2024

It allows reproduction of all figures of the manuscript.


## Figures reproduction code

#### Installation requirements

The visualization code uses Python 3.

A number of standard scientific python packages are needed for the numerical simulations and visualizations. An easy way to install all of these is to install a Python distribution such as [Anaconda](https://www.anaconda.com/):

- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/stable/index.html)


#### Files organization/running the code

All files necessary to reproduce the figures are in the required relative positions.
The python scripts that run the models are in [models](./models), where [lib](./lib) contains some plotting cosmetics definitions, and [plots](./figs_paper_coarse_gr_clean/python_code/plots) contains the scripts to run to produce the figures in the corresponding folder. Pre-processed data for Figs 2,3 and 4 can be found in the subdirectories of [fig2](./figs_paper_coarse_gr_clean/fig2).

For some figures cosmetic changes were done in inkscape as a postprocessing step. In these cases the figures will not be reproduced precisely.

## Simulation and pre-processing code


The simulation and analysis code uses a combination of C, C++, Bash and Python 2.7+. 

The C++ code depends on the packages [FFTW](http://www.fftw.org/) and [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) that need to be installed on the local machine and properly compiled (see [Makefile](./Cpp_code_github/makefile)).


The code is posted here for transparency, it is not meant to be run on local machines as it is. It may not run across platforms as some synthax is Linux specific and is not necessarily POSIX compliant. 
Some snippets also rely on the specific folder structure on the machine. So if the reader wants to try and use this code they should be careful to change the code to reflect their OS and folder structure. 
It should also be noted that for some parameters choices simulations are computationally heavy both in time and memory requirements, and need to be run with the aid of a computing cluster.

#### Files organization

- In  [bash_scripts_handlers_github](./bash_scripts_handlers_github) there are three bash handlers as examples on how we: 1) sweeped through parameters ([grid_parameters_slurm.sh](./bash_scripts_handlers_github/grid_parameters_slurm.sh)), 2) ran the model C++ code handling extinctions and explosions and output directory structure ([arrayjob_coarse_grained_slurm.sh](./bash_scripts_handlers_github/arrayjob_coarse_grained_slurm.sh)), 3) performed analysis and plots on the model outputs ([plot_all_slurm.sh
](./bash_scripts_handlers_github/plot_all_slurm.sh)). These scripts are meant to be run on a  computing cluster mounting a SLURM queuing engine.
- In  [Cpp_code_github](./Cpp_code_github) there is the core C++ code for the model divided in multiple .cpp and .h (header) files. The [input](./Cpp_code_github/input_zuzia.dat) file is an example of input parameters to be passed to the program. The code saves various output files to track the system evolution over time, and to benchmark convolution algorithms speed and approximation errors. The code handling fast convolution using non-homogeneous Fast Fourier Transforms was adapted from open source code documented on [the NFFT webpage](https://www-user.tu-chemnitz.de/~potts/nfft/fastsum.php), the authors of the original code are properly cited in the manuscript.
- In [analysis](./pyscripts_analysis_github/analysis/) there are the python scripts used to perform analysis on the output of the C++ program. [plot_avg_dyn.py](./pyscripts_analysis_github/analysis/plot_avg_dyn.py) collects some average observables time series, the trace in Fig. 3A was generated by this script. [viral_clustering.py](./pyscripts_analysis_github/analysis/viral_clustering.py) clusters and tracks the viral lineages in antigenic space, the trace in Fig. 3B was generated by this script. [persistence_length.py](./pyscripts_analysis_github/analysis/persistence_length.py) infers the persistence time of the lineage trajectories from the output of [viral_clustering.py](./pyscripts_analysis_github/analysis/viral_clustering.py). [rspeciations.py](./pyscripts_analysis_github/analysis/rspeciations.py) computes the number of lineage splits at varius distance thresholds from the output of [viral_clustering.py](./pyscripts_analysis_github/analysis/viral_clustering.py). [summary_features.py](./pyscripts_analysis_github/analysis/summary_features.py) collects summary statistics reported in the manuscript.

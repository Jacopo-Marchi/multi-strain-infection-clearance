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
In  [bash scripts](./scripts_generate_figures) there are the bash handlers used to sweep through parameters space and call the ODE models integrators (in [models](./models_python)), generating synthetic data and organizing the folders in the required relative positions. [plots](./plots) contains the scripts to produce the figures in once synthetic data are generated.
[lib](./lib) contains some plotting cosmetics definitions.


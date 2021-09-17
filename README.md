# 2d-vertex

## About
This repository contains code that can be used in a 2D Vertex Model simulation of tissue mechanics.  This code is similar to the repository [ubuquitous-drosophila](https://github.com/clintondurney/ubuquitous-drosophila), but lacks an ODE model for the biochemistry. Dependencies have been resolved in this version and the code is easier to adapt to other use cases.

For plotting results, please refer to [vertex-plotting](https://github.com/clintondurney/vertex-plotting). 

Functionality is provided as is and written for my own applications.  Little knowledge of Python would be needed to adapt the code to your own use case.  If you are having trouble adapting it, please contact me for help.

## Author:
* [Clinton H. Durney](https://clintondurney.github.io/) (clinton.durney@jic.ac.uk)

## To use:
__Installation__
Use [Conda](https://docs.conda.io/en/latest/) to create the environment from the .yaml file:
```
conda env create -f vertex.yaml python=3.7
```

__Running__
1. Activate the conda environment
```
conda activate vertex
```
2. Run main.py 
```
python main.py
```

This will populate the directory with .pickle files that contain the attributes of the network state and .npy files that save the nodes of the cells circumference oriented counter clockwise. These files can be analyzed using a Jupyter Notebook or visualized using the repository referenced above.

## Citation:
If you find the code useful, please consider citing either or both of the following:
1. [Durney et. al. (2018)](https://www.sciencedirect.com/science/article/pii/S0006349518311615)
2. [Durney et. al. (2021)](https://iopscience.iop.org/article/10.1088/1478-3975/abfa69/meta)









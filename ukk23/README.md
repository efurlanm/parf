# Ukk23

The Notebook `ukk23test01-train.ipynb` generates files containing the neural network (NN) model that is later used in the RTE+RRTMGP-NN model. The implementation uses TensorFlow and Python, and Fortran routines are used to generate the training data set. The idea is to replace the RRTMGP lookup tables with NN.

Based on:

Ukkonen, P., & Hogan, R. J. (2023). Implementation of a machine-learned gas optics parameterization in the ECMWF Integrated Forecasting System: RRTMGP-NN 2.0. Geoscientific Model Development, 16(11), 3241–3261. https://doi.org/10.5194/gmd-16-3241-2023

## CODE AND DATA

RTE+RRTMGP-NN is available on GitHub

- https://github.com/peterukk/rte-rrtmgp-nn (last access: 8 June 2023)
- https://doi.org/10.5281/zenodo.7413935 (Ukkonen, 2022c)
    - https://zenodo.org/records/7413935
    - peterukk/rte-rrtmgp-nn: 2.0, Zenodo [code]

The Fortran programs and Python scripts used for data generation and model training are found in the subdirectory 

- examples/rrtmgp-nn-training

The training data and archived version of RTE+RRTMGP-NN 2.0 with its training scripts can be accessed at 

- https://doi.org/10.5281/zenodo.6576680 (Ukkonen, 2022d)
    - https://zenodo.org/records/7413952
    - Code and extensive data for training neural networks for radiation, used in “Implementation of a machine-learned gas optics parameterization in the ECMWF Integrated Forecasting System: RRTMGP-NN 2.0”, Zenodo [code and data set]

Finally, the optimized version of the ecRad radiation scheme integrated with RRTMGP-NN 2.0 can be accessed at

- https://doi.org/10.5281/zenodo.7148329 (Ukkonen, 2022e)
    - https://zenodo.org/records/7852526
    - Optimized version of the ecRad radiation scheme with new RRTMGP-NN gas optics, Zenodo [code]



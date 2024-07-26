# Ukk23

This sub-repo contains my notes and studies on RTE+RRTMGP-NN.

Mainly based on:

Ukkonen, P., & Hogan, R. J. (2023). Implementation of a machine-learned gas optics parameterization in the ECMWF Integrated Forecasting System: RRTMGP-NN 2.0. Geoscientific Model Development, 16(11), 3241–3261. https://doi.org/10.5194/gmd-16-3241-2023


## NOTEBOOKS

- **`ukk23test01-train.ipynb`**: generates files containing the neural network (NN) model that is later used in the RTE+RRTMGP-NN model. The implementation uses TensorFlow and Python, and Fortran routines are used to generate the training data set. The idea is to replace the RRTMGP lookup tables with NN.

- **`ukk23test01-rfmip-clear-sky.ipynb`**: runs the RFMIP-CLEAR-SKY example.

- **`ecrad01-gprof.ipynb`**: gprof of ecrad executable from ecrad dir.

- **`ukk23eo01-gprof.ipynb`**:  gprof of ecrad executable from ukk23eo01 dir.


## DIRECTORIES

- **`ukk23test01`**: based on <https://zenodo.org/records/7413935>

- **`ukk23eo01`**: based on <https://github.com/peterukk/ecrad-opt>

- **`ecrad`**: based on <https://github.com/ecmwf-ifs/ecrad>


## CODE AND DATA

The code and data come from various sources, such as:

- <https://zenodo.org/records/7413935>
- <https://zenodo.org/records/4030436>
- <https://doi.org/10.5281/zenodo.5833494>
- <https://zenodo.org/records/7413952>
- <https://zenodo.org/records/7852526>

The repositories in <https://github.com/peterukk> also contain code and data, which are distributed across different branches.

RTE+RRTMGP-NN is available on GitHub:

- https://github.com/peterukk/rte-rrtmgp-nn (last access: 8 June 2023)
- https://doi.org/10.5281/zenodo.7413935 (Ukkonen, 2022c) [1]
    - https://zenodo.org/records/7413935 [code]

The Fortran programs and Python scripts used for data generation and model training are found in the subdirectory

- examples/rrtmgp-nn-training

The training data and archived version of RTE+RRTMGP-NN 2.0 with its training scripts can be accessed at

- https://doi.org/10.5281/zenodo.6576680 (Ukkonen, 2022d) [2]
    - https://zenodo.org/records/7413952 [code and data set]

Finally, the optimized version of the ecRad radiation scheme integrated with RRTMGP-NN 2.0 can be accessed at

- https://doi.org/10.5281/zenodo.7148329 (Ukkonen, 2022e) [3]
    - https://zenodo.org/records/7852526 [code]


## NOTES

- The "*.nc" files are from NetCDF4 and their structure can be viewed using the [ToolsUI](https://docs.unidata.ucar.edu/netcdf-java/current/userguide/reading_cdm.html).


## REFERENCES

[1] Ukkonen, P., Pincus, R., Hillman, B. R., Norman, M., fomics, & Heerwaarden, C. van. (2022c). peterukk/rte-rrtmgp-nn: 2.0 (2.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.7413935

[2] Ukkonen, P. (2022). Code and extensive data for training neural networks for radiation, used in “Implementation of a machine-learned gas optics parameterization in the ECMWF Integrated Forecasting System: RRTMGP-NN 2.0”" [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.7413952

[3] Ukkonen, P. (2022). Optimized version of the ecRad radiation scheme with new RRTMGP-NN gas optics [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.7852526

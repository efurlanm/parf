# PINN

[Physics-Informed Neural Networks (PINNs)](https://en.wikipedia.org/wiki/Physics-informed_neural_networks) are a type of universal function approximators that can incorporate the knowledge of physical laws that govern a given set of data in the learning process and that can be described by partial differential equations. This directory contains some experiments I've done with PINN, including references to other repos.

## Files in the current directory

* [burgers-tf1-gpu.ipynb](burgers-tf1-gpu.ipynb) - TensorFlow v1 PINN, adapted from the original [Raissi PINN work](https://github.com/maziarraissi/PINNs/blob/master/appendix/continuous_time_inference%20\(Burgers\)/Burgers.py), running on GPU
* [burgers-tf2-gpu.ipynb](burgers-tf2-gpu.ipynb) - TensorFlow v2 implementation of [Raissi PINN](https://github.com/maziarraissi/PINNs/blob/master/appendix/continuous_time_inference%20\(Burgers\)/Burgers.py), based on the [work of Morzaria](https://github.com/deepmorzaria/Physics-Informed-Neural-Network-PINNs---TF-2.0/blob/master/PINNs_2.ipynb), running on GPU

## Some info available in my other repos

* In [CAP-421 course notes](https://github.com/efurlanm/421) I compared the performance of two implementations of the Burgers viscosity equation, one using a numerical method, and the other using PINN :
  
  * [burgers-tf1_with_graphs.ipynb](https://github.com/efurlanm/421/blob/main/project/burgers-tf1_with_graphs.ipynb) - contains the PINN implementation, based on the [work of Raissi](https://github.com/maziarraissi/PINNs/tree/master/appendix/continuous_time_inference%20(Burgers))
  * [burgers-tf1.ipynb](https://github.com/efurlanm/421/blob/main/project/burgers-tf1.ipynb) - the same implementation, adapted to run on the Santos Dumont supercomputer
  * [burgers-f90.ipynb](https://github.com/efurlanm/421/blob/main/project/burgers-f90.ipynb) - contains the numeric implementation, based on the [work of Bukardt](https://people.sc.fsu.edu/~jburkardt/f_src/burgers_solution/burgers_solution.html)

* In [CAP-418 course notes](https://github.com/efurlanm/418) I made this Notebook :
  
  * [pinn-quadrature-burgers.ipynb](https://github.com/efurlanm/418/blob/main/pinn-quadrature-burgers.ipynb) - based on the works of [Pierre Jacquier](https://github.com/pierremtb/PINNs-TF2.0), [Maziar Raissi](https://github.com/maziarraissi/PINNs), and [John Burkardt](https://people.sc.fsu.edu/~jburkardt/)

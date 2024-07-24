# My personal notes on Machine Learning and other random topics

*Last edited: 2024-04-23*

This repo contains my collection of materials and random notes I take while researching and playing with Machine Learning (ML) and Artificial Neural Networks (NN). It includes topics such as numerical methods, differential calculus, artificial neural networks, libraries, implementations, materials that I ended up using during my research, and other random topics that I also found interesting. It is a work in progress and subject to constant change.


## Contents

Some materials are in the form of simple files and others are organized in subdirectories (unordered list):

* [sparsedynamics](sparsedynamics) - mirror of [sparsedynamics.zip](http://faculty.washington.edu/sbrunton/sparsedynamics.zip) from [Brunton et al.  (2016)](https://doi.org/10.1073/pnas.1517384113), used e.g. in SINDy [[1]](https://rajdandekar.github.io/SINDYExamples_Julia/), [[2]](https://pysindy.readthedocs.io/en/stable/examples/3_original_paper/example.html).
* [nektar](nektar) - [Nektar++](https://www.nektar.info/) Spectral / HP Element Framework running on Colab using [udocker](https://indigo-dc.github.io/udocker/).
* [pinn](pinn) - notes and materials regarding Physics-Informed Neural Networks (PINN).
* [parf](parf) - Random Forest (RF) algorithm in Fortran.
* [rforest](rforest) - RF algorithm in Python.
* [burgers](burgers) - notes related to the convection diffusion equation that is used in several examples in this repo.
* [deepxde](deepxde) - material direct related to the [DeepXDE library](https://deepxde.readthedocs.io).
* [weka](weka) - contains some links about [Weka](https://www.cs.waikato.ac.nz/ml/weka/).
* [horovod](horovod) - (see below).
* [loss](loss) - visualizing the loss landscape of a NN.
* [ukk23](ukk23) - fork of [rte-rrtmgp-nn](https://github.com/peterukk/rte-rrtmgp-nn). The idea is to replace the [RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp) lookup table with NN.


## Horovod

[Horovod](https://horovod.readthedocs.io/en/stable/) was created internally at Uber to make it easy to use a single-GPU training script and successfully scale it to train on many GPUs in parallel.

The horovod directory contains some Notebooks with examples:

* [install-horov-tf1-sd.ipynb](horovod/install-horov-tf1-sd.ipynb) - Installing Horovod and TensorFlow v1 on SDumont supercomputer.
* [hv-tf1-mnist.ipynb](horovod/hv-tf1-mnist.ipynb) - MNIST with TensorFlow v1 using MPI through Horovod, running on SDumont.


## Some info available in my other repos

* In [My MSc repo](https://github.com/efurlanm/msc22) dedicated to my master's thesis, I trained a convolutional NN :

  * The [other](https://github.com/efurlanm/msc22/tree/main/other) directory contains a [PyTorch example](https://github.com/efurlanm/msc22/blob/main/other/pytorch.ipynb) of convolutional NN training using the MNIST database, running on the Santos Dumont supercomputer, [adapted from IDRIS](http://www.idris.fr/eng/jean-zay/gpu/jean-zay-gpu-torch-multi-eng.html).

* In [CAP-351 course notes](https://github.com/efurlanm/351) I made these Notebooks :

  * [project1-mlp.ipynb](https://github.com/efurlanm/351/blob/main/project1-mlp.ipynb) - Multilayer Perceptron (MLP) is a fully connected class of feed-forward artificial neural network (NN).
  * [project2-som.ipynb](https://github.com/efurlanm/351/blob/main/project2-som.ipynb) - a self-organizing map or self-organizing feature map is an unsupervised machine learning technique used to produce a low-dimensional representation of a higher dimensional data set while preserving the topological structure of the data.
  * [project3-vae.ipynb](https://github.com/efurlanm/351/blob/main/project3-vae.ipynb) - in machine learning, a variational auto-encoder, is an artificial neural network architecture introduced by Diederik P. Kingma and Max Welling, belonging to the families of probabilistic graphical models and variational Bayesian methods.
  * [project4-cnn.ipynb](https://github.com/efurlanm/351/blob/main/project4-cnn.ipynb) - a Convolutional Neural Network (CNN, or ConvNet) is a class of artificial neural network (NN), most commonly applied to analyze visual imagery.
  * [project5-rnn.ipynb](https://github.com/efurlanm/351/blob/main/project5-rnn.ipynb) - a Recurrent Neural Network (RNN) is a class of artificial neural networks where connections between nodes can create a cycle, allowing output from some nodes to affect subsequent input to the same nodes.


<br>
<table>
  <tr>
    <td><img src="img/construction.gif"></td>
    <td>This repo is permanently under construction, so its content changes constantly.</td>
  </tr>
</table>

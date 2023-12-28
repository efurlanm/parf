# My personal notes on machine learning

*Last edited: 2023-10-15*

This repository contains my collection of materials and random notes that I take while researching and playing Machine Learning (ML). It is a work in progress and is subject to constant change.

## Contents of this repo

Materials on the various subjects are organized into sub-directories :

* [pinn](pinn) -  notes and materials regarding Physics-Informed Neural Networks (PINN)
* [parf](parf) - Random Forest (RF) algorithm in Fortran
* [rforest](rforest) - RF algorithm in Python
* [burgers](burgers) - notes related to the convection diffusion equation that is used in several examples in this repo
* [deepxde](deepxde) - material direct related to the [DeepXDE library](https://deepxde.readthedocs.io)
* [datasets](datasets) - some datasets used in the works

## Horovod

[Horovod](https://horovod.readthedocs.io/en/stable/) was created internally at Uber, as a model for distributed use of TensorFlow. The horovod directory contains some Notebooks with examples.

* [hv-tf1-mnist.ipynb](horovod/hv-tf1-mnist.ipynb)

## Some info available in my other repos

* In [My MSc repo](https://github.com/efurlanm/msc22) dedicated to my master's thesis, I trained a convolutional ANN :
  
  * The [other](https://github.com/efurlanm/msc22/tree/main/other) directory contains a [PyTorch example](https://github.com/efurlanm/msc22/blob/main/other/pytorch.ipynb) of convolutional ANN training using the MNIST database, running on the Santos Dumont supercomputer, [adapted from IDRIS](http://www.idris.fr/eng/jean-zay/gpu/jean-zay-gpu-torch-multi-eng.html)

* In [CAP-351 course notes](https://github.com/efurlanm/351) I made these Notebooks :
  
  * [project1-mlp.ipynb](https://github.com/efurlanm/351/blob/main/project1-mlp.ipynb) - Multilayer Perceptron (MLP) is a fully connected class of feed-forward artificial neural network (ANN)
  * [project2-som.ipynb](https://github.com/efurlanm/351/blob/main/project2-som.ipynb) - a self-organizing map or self-organizing feature map is an unsupervised machine learning technique used to produce a low-dimensional representation of a higher dimensional data set while preserving the topological structure of the data
  * [project3-vae.ipynb](https://github.com/efurlanm/351/blob/main/project3-vae.ipynb) - in machine learning, a variational auto-encoder, is an artificial neural network architecture introduced by Diederik P. Kingma and Max Welling, belonging to the families of probabilistic graphical models and variational Bayesian methods
  * [project4-cnn.ipynb](https://github.com/efurlanm/351/blob/main/project4-cnn.ipynb) - a Convolutional Neural Network (CNN, or ConvNet) is a class of artificial neural network (ANN), most commonly applied to analyze visual imagery
  * [project5-rnn.ipynb](https://github.com/efurlanm/351/blob/main/project5-rnn.ipynb) - a Recurrent Neural Network (RNN) is a class of artificial neural networks where connections between nodes can create a cycle, allowing output from some nodes to affect subsequent input to the same nodes

## Random notes

* I use [Zotero](https://www.zotero.org/) as a research assistant, and the library I'm using is in the file [ml.bib](ml.bib) in BibTeX format.

<br>
<table>
  <tr>
    <td><img src="img/construction.gif"></td>
    <td>This work is permanently under construction, so its content changes constantly.</td>
  </tr>
</table>

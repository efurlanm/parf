# Poisson TF PT PD

TF1 = TensorFlow v1

TF2 = TensorFlow v2

PT = PyTorch
PD = Paddle

Comparison of different backends, TensorFlow, PyTorch, and Paddle, for the  [Poisson equation in L-shaped domain](https://deepxde.readthedocs.io/en/latest/demos/pinn_forward/poisson.lshape.html), running on a GTX 1050 3 MB GPU:

- [poisson-tf1-01.ipynb](poisson-pt-01.ipynb) - using TensorFlow v1 (TF1) backend
- [poisson-tf2-01.ipynb](poisson-pt-01.ipynb) - using TensorFlow v2 (TF2) backend
- [poisson-pt-01.ipynb](poisson-pt-01.ipynb) - using PyTorch (PT) backend
- [poisson-pd-01.ipynb](poisson-pd-01.ipynb) - using Paddle (PD) backend

Training times, first using the Adam optimizer, then retraining using L-BFGS to get a smaller loss. The best performing backend is TensorFlow v1:

| backend | Adam (s) | L-BFGS (s) |
|:-------:|:--------:|:----------:|
| TF1     | 167      | 53         |
| TF2     | 209      | 754        |
| PT      | 395      | 311        |
| PD      | 344      | 409        |


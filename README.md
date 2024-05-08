# VAMP with Arbitrary IID Noise Priors
##### Mohamed Akrout, Tiancheng Gao, Faouzi Bellili, Amine Mezghani
This repository contains the Matlab code of the algorithm proposed in our ICASSP 2024's paper: [Vector Approximate Message Passing With Arbitrary I.I.D. Noise Priors](https://arxiv.org/abs/2402.04111).


## Abstract
Approximate message passing (AMP) algorithms are devised under the Gaussianity assumption of the measurement noise vector. In this work, we relax this assumption within the vector AMP (VAMP) framework to arbitrary independent and identically distributed (i.i.d.) noise priors. We do so by rederiving the linear minimum mean square error (LMMSE) to accommodate both the noise and signal estimations within the message passing steps of VAMP. Numerical results demonstrate how our proposed algorithm handles non-Gaussian noise models as compared to VAMP. This extension to general noise priors enables the use of AMP algorithms in a wider range of engineering applications where non-Gaussian noise models are more appropriate.

## Repository Structure
This repository contains two folders:
  - **priors**: it contains the code of the supported priors: Bernoulli-Gaussian (bg) and binary priors.
  - **lmmse**: it contains the code of the LMMSE estimation steps of the signal $\boldsymbol{x}$ and the noise $\boldsymbol{w}$.

The repository also contains two files:
  - **VampNoiseIID.m**: it contains the VAMP algorithm supporting i.i.d. noise priors which uses the functions in the folders `priors` and `lmmse`.
  - **main.m**: it contains the code of the main file that runs the algorithm.

## Citing the paper (bib)

If you make use of our code, please make sure to cite our paper:
```
@INPROCEEDINGS{10446747,
  author={Akrout, Mohamed and Gao, Tiancheng and Bellili, Faouzi and Mezghani, Amine},
  booktitle={ICASSP 2024 - 2024 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)}, 
  title={Vector Approximate message Passing with Arbitrary I.I.D. Noise Priors}, 
  year={2024},
  pages={9596-9600},
  doi={10.1109/ICASSP48485.2024.10446747}
}
```

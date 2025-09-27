# TSP_PSU_CCPD

Power Spectrum Unmixing via Coupled Tensor Decomposition

This repository contains the code used in the paper:

Chen, Xiaonan and Wang, Jun. "Power Spectrum Unmixing via Coupled Tensor Decomposition: Unified Identifiability and Algebraic Solution." IEEE Transactions on Signal Processing, 2025.

Citation

```
@article{chen2025power,
  title={Power Spectrum Unmixing via Coupled Tensor Decomposition: Unified Identifiability and Algebraic Solution},
  author={Chen, Xiaonan and Wang, Jun},
  journal={IEEE Transactions on Signal Processing},
  year={2025},
  publisher={IEEE}
}
```

Contents
- MATLAB functions and scripts used in the experiments
- `functions/` — algorithm implementations and supporting code
- `verifications_of_theorem/` — experimental scripts for the paper's theorems
- `adding noise/` — experimental scripts under noisy settings. We use two kind of denoising approach. Functions with suffix "nd" refer to the non-diagonal sub-tensors denoising, see more details in Appendix D

SDR data is provided by [Prof. Xiao Fu](https://web.engr.oregonstate.edu/~fuxia/index.html).
[ref](https://ieeexplore.ieee.org/document/7175044) X. Fu, N. D. Sidiropoulos, J. H. Tranter, and W.-K. Ma, “A factor
analysis framework for power spectra separation and multiple emitter
localization,” IEEE Trans. Signal Process., vol. 63, no. 24, pp. 6581–
6594, 2015

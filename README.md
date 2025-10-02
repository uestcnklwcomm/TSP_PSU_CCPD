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
- `adding noise/` — experimental scripts under noisy settings. Scripts with suffix "_nd_" refer to the implementation of non-diagonal sub-tensors denoising technique, see more details in Appendix D
- `sdr_experiments/` — Scripts of laboratory experiments using software defined radios. The data is provided by [Prof. Xiao Fu, OSU](https://web.engr.oregonstate.edu/~fuxia/index.html), [ref](https://ieeexplore.ieee.org/document/7175044), [download link](https://drive.google.com/drive/folders/1wtxkoNKCkIHH8Lm9BSjH7CvCPnEXr7rP?usp=drive_link)
- `Baselines/`
- `Case 1)/` — (1) successive projection algorithm (SPA), [ref](https://ieeexplore.ieee.org/document/7463032) X. Fu, N. D. Sidiropoulos, and W.-K. Ma, “Power spectra separation via
structured matrix factorization,” IEEE Trans. Signal Process., vol. 64,
no. 17, pp. 4592–4605, Sep. 2016.   (2) ESRPIT [ref](https://ieeexplore.ieee.org/document/7175044) X. Fu, N. D. Sidiropoulos, J. H. Tranter, and W.-K. Ma, “A factor
analysis framework for power spectra separation and multiple emitter
localization,” IEEE Trans. Signal Process., vol. 63, no. 24, pp. 6581–
6594, 2015.

- `Case 2)/` — Carrier-frequency-based Compressed sensing, [ref](https://ieeexplore.ieee.org/document/9674222) H. Wang, J. Fang, H. Duan, and H. Li, “Compressive wideband spectrum
sensing and signal recovery with unknown multipath channels,” IEEE
Trans. Wirel. Commun., vol. 21, no. 7, pp. 5305–5316, Jul. 2022.
  

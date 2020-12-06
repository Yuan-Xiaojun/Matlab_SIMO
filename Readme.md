# Passive Beamforming and Information Transfer via Large Intelligent Surface

This package contains the official implementation of the **passive beamforming method** and the **two-step signal detection methods** proposed in the paper: 

> W. Yan, X. Yuan and X. Kuai, "Passive Beamforming and Information Transfer via Large Intelligent Surface," *IEEE Wireless Communications Letters*, vol. 9, no. 4, pp. 533-537, Apr. 2020, doi: [10.1109/LWC.2019.2961670] (https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8941126).

## Introduction

This paper proposed a novel **passive beamforming and information transfer (PBIT)** technique by adopting spatial modulation on the index of the LIS (or RIS) elements, in which the LIS simultaneously enhances the primary communication (by passive beamforming) and sends its private data to the receiver (by spatial modulation). In this paper, we developed a passive beamforming method to improve the average receive signal-to-noise ratio (SNR), and established a two-step approach at the receiver to retrieve the information from both the transmitter and the LIS. 


## Code Structure

`main.m`: Set system parameters and output simulation results. 

`SIMO_algorithm`: Generate system model.

`CVX_Optimal`: Optimize RIS phase shift.

`Constell_Modulate` and `Constell_Mapping`: Constellation modulation and demodulation.

`BiGAMP`: BiGAMP algorithm to detect user signals.

`GAMP_Bernoulli`: GAMP algorithm to detect RIS data.

`SIMO_result.txt`: Save the simulation results. 

## Citation
```
@article{yan2019passive,
  title={Passive beamforming and information transfer via large intelligent surface},
  author={Yan, Wenjing and Yuan, Xiaojun and Kuai, Xiaoyan},
  journal={IEEE Wireless Communications Letters},
  volume={9},
  number={4},
  pages={533--537},
  month={Apr.}
  year={2020},
  publisher={IEEE}
}
```





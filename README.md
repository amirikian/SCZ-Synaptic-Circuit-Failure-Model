# SCZ-Synaptic-Circuit-Failure-Model

This repository provides the codes for the spiking network model simulation and the mean field and linear stability analysis presented in the paper: [A prefrontal network model operating near steady and oscillatory states links spike desynchronization and synaptic deficits in schizophrenia](), Crowe DA, Willow A, Blackman RK, DeNicola AL, Chafee MV, Amirikian B, *eLife* 2023. 

The repository contains three directories:

1. StateDiagrams
2. SpikingNtwrks
3. SpikeTrains

The codes are implemented in MATLAB R2018b.

## 1. StateDiagrams
This directory contains the code that solves the mean field and linear stability equations and computes network state diagrams. The code reproduces Fig.2 A, B, the state diagrams shown in Fig.5 A, B, and saves computed synaptic conductance parameters for the steady and critical state primary networks in the `Parameters` subdirectory.

## 2. SpikingNetworks
This directory contains the code that simulates the dynamics of recurrent spiking network model for given synaptic conductance parameters. Using the conductance parameters for the steady and critical state primary networks from the `Parameters` subdirectory, the code generates spike trains of individual neurons corresponding to the networks in Fig.3 and saves them in the `Simulations` subdirectory. 

## 3. SpikeTrains
This directory contains the code that, using the spike trains from the `Simulations` subdirectory, reproduces panels A1, B1, A2, B2 in Fig.3, and plots the corresponding spike correlation functions as those in Fig.4.
## Overview
In this discussion we will look at a non-reversible Monte Carlo simulation. Most Monte Carlo algorithms require the transition probabilities to satisfy the detailed balance equation. But especially for Markov chains on finite state space we can easily proof that we only need irreducibility and aperiodicity to achieve convergence to a unique invariant distribution. That is why we will relax the detailed balance condition and construct an algorithm with a non-reversible state update to speed up convergence.

![Screenshot Example](Matlab%20Code/Example%20Plots/nonrevMFI1D.png)

 
## Information
Please refer to the documents in the folder 'Materials' for further information. A quick introduction can be found in the document 'Handout'. For a more visual introduction consider to read the document 'Slides'. The document 'Convergence Theorem Proof' can be seen as extra material.

## Setup and Execution
This repository contains three examples which can be executed with the following parameters. Warning: execution time might be long.

**Example (Non-reversible 1D-Ising-Model over M)**

nonrevMFI1D(theta,plotRun,calcIndicators,samples)

``` nonrevMFI1D(0.001,1,0,10000) ```

**Example (Non-reversible 2D-Ising-Model over M)**

nonrevMFI2DM(theta,plotRun,calcIndicators,samples)

``` nonrevMFI2DM(0.001,1,0,10000) ```

**Example (Non-reversible 2D-Ising-Model over E)**

nonrevNNI2DE(theta,plotRun,calcIndicators,samples)

``` nonrevNNI2DE(0.001,1,0,10000) ```

## Feedback and Support
Please feel free to open an issue if you find a bug or seek support. This project probably will not have any further development. That is why I encourage you to fork the project.

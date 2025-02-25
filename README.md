# ND-joint-distribution-simulator
## Overview
The ND Joint Distirbution Simulator creates a random joint distribution of any length and dimensionality (within reason) based on an input joint distribution, such as one obtained from observations.  It was developed to support Monte Carlo modeling experiments of atmospheric conditions, where there are likely relationships between the variables.  
## Installation
This Matlab package consists of a wrapper script which prepares the input data and produces diagnostic plots (`sampleND_wrapper.m`) and the distribution generation function itself (`inverseTransformSampleND.m`).  Copy the files into the desired place in your matlab path.  Matlab's `Deep Learning Toolbox` must also be installed.
## Usage
The wrapper script includes 3-dimensional synthetic data for illustration.  The three data vectors are then put into a `var` stucture for downstream processing.  The user will change this `var` structure to contain their desired data (of whatever dimension).  Each vector of the `var` structure must be of the same length.

There are 2 input parameters that are set inside the wrapper script.
- The output data samples are produced on a discrete basis.  The vector `nBins` controls how many bins to divide each input data vector into.  These bins are equally-spaced between the min and max of each data vector. The length of `nBins` must be the same as the number of dimensions in the input data.
- The number of desired joint output samples `nSamples`.
The wrapper script then computes bin edges, the centered values to be associated with each bin, and an N-dimensional probability matrix based on the input data.  These and the input parameters are passed to the `inverseTransformSampleND` function.

The `inverseTransformSampleND` function computes a cumulative distribution function based on the input probabilty matrix and applies Inverse Transform Sampling to marginal distributions along each dimension to extract random values.  These marginal distributions are computed across all dimensions up to the one being processed, and their expression is hardcoded for up to 10-dimensional input for efficiency.  For higher dimensions the expression is dynamically constructed.  The function then returns a sample structure `samp` of the desired length with elements corresponding to the input structure `var`.

The wrapper then turns the output structure `samp` to table `sampT`.  The external writing of output is left to the user.  The wrapper also produces diagnostic plots of input vs sample marginal 1D histograms as well input vs sample 2D histograms for all possible 2-way variable combinations.

This code has been employed to produce up to 10M samples of 5-dimensional joint distributions using 100 bins across 4 of the dimensions and 50 bins for the fifth.  The primary constraint is the memory size of the N-dimensional probability distribution (100<sup>4</sup> * 50 = 37.2 GB for the aforementional case).  The user will need to trade off between the number of bins and number of dimensions to satisfy their own memory constraints.

Creating 10M samples for the above case took approximately 20 minutes of wall time on an Intel Core i9-11950H processor (8 cores) running at a base clock speed of 2.60 GHz.  No attempt has yet been made to parallelize the code.
## Disclaimer
DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
 
This material is based upon work supported by the Federal Aviation Administration under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the Federal Aviation Administration.

© 2024 Massachusetts Institute of Technology.
 
Subject to FAR52.227-11 Patent Rights - Ownership by the contractor (May 2014)

The software/firmware is provided to you on an As-Is basis

Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb2014). Notwithstanding any copyrightnotice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.

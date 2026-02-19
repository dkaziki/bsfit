The bsfit toolbox provides a principled framework for reducing and interpreting high-dimensional cross-bispectral EEG data. By using a low-rank tensor decomposition, the toolbox expresses complex sensor-level interactions as a product of a single spatial mixing matrix and a compact source-interaction tensor.
​
Directory of Functions:

Core Algorithm and Data Processing

• bsfit.m: The primary function for fitting the low-rank model to a bispectral tensor using the Levenberg-Marquardt (LM) algorithm.
• data2source.m: Projects sensor-level data into source space to identify the spatial origin of bispectral components.
• data2bs_event.m: Specifically designed for computing cross-bispectral tensors for event-related EEG data.
• data_single.mat: A sample dataset provided for testing and validation.

Example Pipelines

• main_bsfit_short.m: A simplified execution script for rapid validation of the fitting process.
• main_bsfit_long.m: A comprehensive end-to-end pipeline covering the full analysis from data to source-level interpretation.
• run_sim_example.m: A template for generating realistic multichannel EEG data

Installation & Requirements
• Language: MATLAB 
• Dependencies: Requires the METH toolbox for forward/inverse calculations, visualization, and MOCA implementation.

Citation
If you use this toolbox, please cite:
Kaziki, D., Engel, A. K., & Nolte, G. (2025). Low-Rank Tensor Decomposition for Cross-Bispectral Analysis of EEG Data. bioRxiv.

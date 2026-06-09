# HamNet-Ensemble

## Overview

This repository contains the code used in the paper *[XXX]* (DOI: XXX), designed to help reproduce the research. The project implements a Deep Ensemble method based on the [DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3.git) framework. The ensemble mechanism is directly written into the original DeepH-E3 code without any modifications to the graph neural network architecture.

> **Note**: The modified sections in `kernel.py` and `parse_configs.py` are marked with comments `# DE-DeepH: ...` for easy identification.

## Based On

This codebase implements Deep Ensemble methods for Hamiltonian prediction, using [DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3.git) as an experimental framework for implementation reference.

## Project Structure

The repository includes the following files and directories:

- `kernel.py`: Core ensemble module (modified from original DeepH-E3)
- `parse_configs.py`: Configuration file parser (modified from original DeepH-E3)
- `Bilayer_graphene_eval_ensemble.ini`: Reference evaluation configuration
- `tools/`: Utility scripts directory containing tools for generating original DFT inputs, post-processing prediction data, and plotting results

### Core File Descriptions

- **kernel.py**: Implements the core logic for ensemble evaluation, directly embedded into the DeepH-E3 inference workflow. All modifications are annotated with `# DE-DeepH:` comments.
- **parse_configs.py**: Extended configuration parser supporting ensemble-related parameter settings. All modifications are annotated with `# DE-DeepH:` comments.
- **Bilayer_graphene_eval_ensemble.ini**: An example evaluation configuration file for bilayer graphene, ready to use as a template.
- **tools/**: Contains utility scripts for generating DFT inputs, post-processing predictions, and plotting results.

## Installation and Usage

### Requirements

This project is built upon [DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3.git). Please ensure it is correctly installed.

### Integration Steps

1. **Replace Core Files**
   Copy `kernel.py` and `parse_configs.py` from this project to your DeepH-E3 installation directory, overwriting the original files:

       cp kernel.py path/to/your/DeepH-E3/deephe3/kernel.py
       cp parse_configs.py path/to/your/DeepH-E3/deephe3/parse_configs.py

2. **Configure Evaluation Parameters**
   Refer to the provided `Bilayer_graphene_eval_ensemble.ini` file and adjust the evaluation parameters according to your specific needs.

3. **Run Evaluation**
   Start the evaluation task following the standard DeepH-E3 workflow.

## Important Note

**Complex tensor operations are currently not supported.** Scenarios involving complex tensors, such as Hamiltonian matrices with spin-orbit coupling (SOC) effects, are beyond the scope of this work.
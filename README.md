# Beacon Calculus Models for RNA Polymerase II Phosphorylation Dynamics

This repository contains stochastic process algebra models designed to investigate the competing hypotheses of *cis*- vs. *trans*-phosphorylation of RNA Polymerase II (Pol II) in response to DNA damage. The models are developed using the Beacon Calculus process algebra and the `bcs` simulator, and they produce data for comparative analysis against experimental observations.

## Project Overview

In this project, we model the spatiotemporal dynamics of Pol II as it elongates along a gene. The core objective is to understand the mechanism by which Pol II becomes phosphorylated at Serine 7 (Ser7-P) following DNA damage. We hypothesize two distinct mechanisms:

1.  **Cis-Phosphorylation (Baseline Model):** Phosphorylation occurs directly on Pol II molecules that have stalled at a DNA damage site.
2.  **Trans-Phosphorylation (Experimental Model):** A signal is propagated from the site of damage, causing phosphorylation to occur on Pol II molecules at a distinct, distant location.

The models simulate Pol II movement, damage encounter, signal propagation via beacons, and phosphorylation, allowing for a direct comparison of the predicted phosphorylation profiles for each mechanism.

## Models

The repository contains the following Beacon Calculus model files:

* **`cis_model.bc`**: A model where Pol II phosphorylation occurs only at the immediate site of DNA damage. This serves as the baseline for comparison.
* **`trans_model.bc`**: An advanced model that incorporates a `phos` beacon-based signaling mechanism. When a Pol II encounters DNA damage, it launches a `phos` beacon. Another Pol II, at a specific downstream location (`i=2` or 200 bp), can detect this signal and become phosphorylated. This model also uses a `d` flag to ensure that each Pol II instance only emits a damage signal once per encounter, preventing redundant signaling.

## Requirements

To run these simulations and analyze the output, you will need:

* **Beacon Calculus Simulator (`bcs`)**: The `bcs` simulator is required to execute the model files. It can be found at: [https://github.com/MBoemo/bcs.git](https://github.com/MBoemo/bcs.git)
* **Python 3**: For running the data processing and plotting scripts.
* **Python Libraries**: The analysis scripts depend on `numpy`, `matplotlib`, and `scipy`.

## Getting Started

1.  **Download and install the Beacon Calculus Simulator (`bcs`):**
    First, you need to download and compile the `bcs` simulator from the official repository:
    ```bash
    git clone https://github.com/MBoemo/bcs.git
    cd bcs
    git checkout v1.0.1
    make
    ```
    This will create the `bcs` executable in the `bcs/bin` directory. Follow the complete installation instructions at [https://github.com/MBoemo/bcs](https://github.com/MBoemo/bcs) for detailed setup steps.

2.  **Clone this repository:**
    ```bash
    git clone https://github.com/shutongye/bcs_Pol_II_phosphorylation.git
    cd bcs_Pol_II_phosphorylation
    ```
3.  **Run a simulation:**
    Navigate to either the `cis` or `trans` directory and run a simulation. The `-s` flag sets the number of simulations,  `-t` allocates the number of threads to use and `-o` navigates to the output directory.
    ```bash
    # Example for cis model
    /path/to/bcs/bin/bcs -s 1 -t 48 -o my_cis_output.bcs cis_model.bc

    # Example for trans model
    /path/to/bcs/bin/bcs -s 1 -t 48 -o my_trans_output.bcs trans_model.bc
    ```
4.  **Analyze the output:**
    Use the provided Python scripts (e.g., `analysis.py`) to process the `.bcs` output files. These scripts will extract Pol II positions, calculate densities, and generate plots comparing the models.

## License

This project is licensed under the MIT License. This is a permissive license that allows for free use, modification, and distribution of the software, and is a standard for open-source projects.
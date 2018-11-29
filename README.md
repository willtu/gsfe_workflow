# A high-throughput computation framework for generalized stacking fault energies of pure metals
This is the workflow for calculating GSFE (Genralized Stacking Fault Energy) which is based on Materials Project's [pymatgen](https://github.com/materialsproject/pymatgen), [fireworks](https://github.com/materialsproject/fireworks), [custodian](https://github.com/materialsproject/custodian) and [atomate](https://github.com/hackingmaterials/atomate).

## Getting Started
Before using these workflow scripts, you need to install atomate and deploy the cluster, which can refer to the [official document](https://atomate.org/installation.html). 

After installing atomate and  deploying the cluster, you can use this workflow which is made up of three scripts:

* slip_systems.yaml
This script contains typical slip systems for BCC, FCC and HCP structures and one example element for each used to identify the input structure.

* image_generation.py
This script aims to construct the needed structures accoding to input parameters.

* workflow4gsfe.py
This script is used to modify the input and submit the calculation jobs.

Therefore, you just need to modify the workflow4gsfe.py to submit desired GSFE computation. You need to assign the layer (layer of slab structure), mp_id (which could be mp_id from [material projects](https://www.materialsproject.org/) or the local file representing bulk structure), kpoint (k mesh), and other incar settings (including IMAGES, EDIFF, EDIFFG and so on). After completing above settings, you can execute this python script and generate the computation jobs.

## Versioning
This is a beta version, we will update and synchronize it in the future, welcome to follow and fork us, thanks for your attention.

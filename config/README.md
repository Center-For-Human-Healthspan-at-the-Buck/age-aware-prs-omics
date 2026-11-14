# MOVE configuration files

This directory contains the study-specific configuration files used for the
perturbable variational autoencoder (VAE) analyses reported in this study.

The VAE analyses were performed using the publicly available MOVE version 1.5.0:

https://github.com/RasmussenLab/MOVE


Reference:

Allesøe RL et al.
Discovery of drug–omics associations in type 2 diabetes with generative
deep-learning models.
Nature Biotechnology. 2023.
https://doi.org/10.1038/s41587-022-01520-x

## Files

### vae_resid_run__id_assoc_bayes.yaml

Configuration used for the VAE in which age-related variation was removed
from the omic features before model training.

### vae_age__run__id_assoc_bayes.yaml

Configuration used for the VAE in which chronological age and age-stratified
PRS features were included during model training.

The YAML files record the model settings. Age representation was defined during
upstream input preparation and is not encoded in the YAML files.

## Reproducing the study

The R scripts do not retrain the VAE. They use the aggregate association outputs
in `data/vae/` to generate the Main, Supplementary and Extended Data panels.

Researchers wishing to rerun the VAE from participant-level data should install
MOVE from the official repository above and use the corresponding YAML
configuration together with appropriately formatted input data.

Participant-level Arivale data are not distributed because of participant
consent and data-use restrictions.

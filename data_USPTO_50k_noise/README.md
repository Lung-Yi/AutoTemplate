# Generate noisy dataset.
Use the following command to generate the noisy USPTO-50k datasets.

There are 3 types of noise we include in the noisy datasets:
(1) removing one reactant
(2) randomly changing the product atom 
(3) randomly changing the atom mapping 
```
cd scripts
python prepare_noisy_data.py
```
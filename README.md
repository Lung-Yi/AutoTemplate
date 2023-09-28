# AutoTemplate: automatically reaction data curation using reaction template
This is a data-preprocess tool that curates the reaction SMILES.

The manuscript of this repository is in preparation.

![alt text](docs/abstract_TOC.svg)

The experimental results for selected Reaxys dataset is:

![alt text](docs/output.svg)

## Installation
1. `git clone https://github.com/Lung-Yi/AutoTemplate.git`
2. `conda create --name autotemplate python=3.7`
3. `pip install rxnmapper`
4. `pip install rdchiral`
5. `pip install CGRTools`
6. `conda install -c anaconda networkx`
7. `conda install -c anaconda pandas`
8. `pip install openpyxl`
9. `conda install -c conda-forge matplotlib`
10. `pip install py-mini-racer`


## Preprocess the Reaxys dataset
### 1. Prepare the dataset
Check all the directories in ./data_reaxys/

All types of the reaction have their corresponding Reaction ID records in the (.txt) files. Please download the files (.xlsx) on https://www.reaxys.com/#/search/quick

### 2. Preprocess the dataset
(1) use RXNMapper for atom-mapping
set the RXN variable for preprocess reaction: `RXN=AdamsDecarboxylation`
```
python rxnmapper_mapping.py --input_dir data_reaxys/${RXN} \
    --output_file data/${RXN}/MappingResult_${RXN}.txt
```

(2) extract super general reaction templates and apply them on the original reaction
```
python auto_process.py \
    --input_file data/${RXN}/MappingResult_${RXN}.txt \
    --radius 0 \
    --threshold 5
```
## Examine the preprocessed results.
```
python post_analysis.py
```

## Universal reaction template extraction tutorial
For further details about the reaction template extraction, please refer to:

scripts/examples.ipynb
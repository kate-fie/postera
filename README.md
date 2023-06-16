# Postera API Access

API Documentation: https://api.postera.ai/api/v1/docs/#section/Examples

Make sure to set your API Key `"MANIFOLD_API_KEY"`

### Run a SA calc:
Usage: `postera_sa.py [path_to_csv] [index_of_column_of_SMILES_in_csv]`

Returns: list of SA scores for each molecule in csv where 0=difficult to synthesize and 1=readily available 

*NOTE* The SA scores returned are the complement from the Manifold SA scores, where 0=readily available and 1=difficult to synthesize.   

Example:
```
python postera_sa.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1
```

### Run a retrosynthesis calc:
Usage: `postera_retrosynthesis.py [path_to_csv] [index_of_column_of_SMILES_in_csv] [calc_on_first_n]`

Returns: dictionary of retrosynthesis route info for each molecule. i.e. reactants, maxleadtime for reactant, reaction name, etc.
```
python postera_retrosynthesis.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
```

### Run a similarity search:
Usage: `postera_similaritySearch.py [path_to_csv] [index_of_column_of_SMILES_in_csv] [calc_on_first_n]`

Returns: found hits for each query molecule
```
python postera_similaritySearch.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
```

### Run an exact search:
Usage: `postera_exactSearch.py [path_to_csv] [index_of_column_of_SMILES_in_csv] [calc_on_first_n]`

Returns: found hits for each query molecule
```
python postera_exactSearch.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
```
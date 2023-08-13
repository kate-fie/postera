import sqlitedict
import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
import sys


def calculate_tanimoto_similarity(smiles1, smiles2):
    # Create RDKit molecules from SMILES strings
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Generate fingerprints for molecules
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)

    # Calculate Tanimoto similarity
    similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
    #Whats the distribution of the similarity scores?
    print(similarity)

    return similarity

def all_reactions_seen(dictionary, rxn_list):
    rxn_order = []
    price_list = []
    for route in dictionary['routes']:
        rxn_order = [reaction['name'] for reaction in route['reactions']]
        for molecule in route['molecules']:
            price_list.extend([entry.get('bbPriceRange') for entry in molecule['catalogEntries'] if molecule['isBuildingBlock']])
            price_list.extend([entry.get('scrPriceRange') for entry in molecule['catalogEntries'] if not molecule['isBuildingBlock']])

    rxn_list.extend(rxn_order)
    price_list.extend(price_list)
    return rxn_list, price_list

def analyse_routes(dictionary):
    data = []
    scaffold_smiles = 'O=C(C1CN([*:0])C(c2ccc(Cl)cc21)=O)Nc3c(cccc4)c4cnc3'
    encoded_reactions = [
    "Amidation",
    "Amide schotten - baumann",
    "Reductive amination",
    "N-nucleophilic aromatic substitution",
    "Sp2-sp2 Suzuki coupling",
    "Sulfonamide schotten-baumann",
    "Boc protection",
    "Boc deprotection",
    "Buchwald hartwig amination - NMP",
    "Buchwald hartwig amination - EtOH",
    "Buchwald hartwig thiolation",
    "Buchwald hartwig etherification - EtOH",
    "Buchwald hartwig etherification - NMP",
    "Sonogashira coupling - EtOH",
    "Sonogashira coupling - NMP",
    "Heck coupling",
]
    postera_terminology_reactions = [
        'Amide Schotten-Baumann with amine',
        'Buchwald-Hartwig amidation with amide-like nucleophile',
        'Ester amidation',
        'Suzuki coupling',
        'Amidation',
        'Sulfonic ester Schotten-Baumann'
        'Buchwald-Hartwig amidation',
    ]

    for route in dictionary['routes']:
        similar_metric = 0.9
        smiles = [molecule['smiles'] for molecule in route['molecules']]
        # for i, smiles in enumerate(smiles):
        #     if calculate_tanimoto_similarity(smiles, scaffold_smiles) > similar_metric:
        #         similar_smiles.append(smiles)
        #         print(smiles)
        num_steps = len(route['reactions'])
        rxn_order = [reaction['name'] for reaction in route['reactions']]
        reactants = []
        building_blocks = []
        catalog_names = []
        lead_time = []
        similar_smiles = []


        for reaction in route['reactions']:
            reactant_smiles = reaction['reactantSmiles']
            # for i, smiles in enumerate(reactant_smiles):
            #     if calculate_tanimoto_similarity(smiles, scaffold_smiles) > similar_metric:
            #         similar_smiles.append(reactant_smiles)
            reactants.append(tuple(reactant_smiles))
            #catalog_names.extend([entry.get('catalogName') for entry in route['molecules'] if entry['smiles'] in reactant_smiles])

        for molecule in route['molecules']:
            building_blocks.append(molecule['smiles'] if molecule['isBuildingBlock'] else None)
            catalog_names.extend([entry.get('catalogName') for entry in molecule['catalogEntries'] if molecule['isBuildingBlock']])
            lead_time.extend([entry.get('bbLeadTimeWeeks') for entry in molecule['catalogEntries'] if molecule['isBuildingBlock']])
        #building_blocks_for_rxn = [building_blocks[i] if i < len(building_blocks) else None for i in range(num_steps)]
        #catalog_name_for_building_block = [catalog_names[i] if i < len(catalog_names) else None for i in range(num_steps)]

        car_route = all(name in postera_terminology_reactions for name in encoded_reactions)
        #car_route = any(any(name in term for term in postera_terminology_reactions) for name in rxn_order)

        data.append({
            'SMILES': smiles[0],
            'CAR_route': car_route,
            'num_steps': num_steps,
            'rxn_order_first_to_last': rxn_order,
            'reactants': reactants,
            'BuildingBlocks': building_blocks,
            'catalogName_for_BuildingBlock': catalog_names,
            #'max_lead_time': max(lead_time)
            #f'similar_smiles_{similar_metric}': similar_smiles
        })

    return pd.DataFrame(data)

# Only looking at 1 compound
one_compound = True
i = 0
rxn_list = []
df = pd.DataFrame()
similarity_df = pd.DataFrame()
for key, value in table.items():
    i += 1
    dictionary = value
    if i == 1:
        df = analyse_routes(dictionary)
    else:
        df = df.append(analyse_routes(dictionary))
    rxn_list, price_list = all_reactions_seen(dictionary, rxn_list)

#     # Calculate and store Tanimoto similarity values
#     smiles = [molecule['smiles'] for route in dictionary['routes'] for molecule in route['molecules']]
#     scaffold_smiles = 'O=C(C1CN([*:0])C(c2ccc(Cl)cc21)=O)Nc3c(cccc4)c4cnc3'
#     similarity_values = [calculate_tanimoto_similarity(smiles, scaffold_smiles) for smiles in smiles]
#     similarity_df = similarity_df.append(pd.DataFrame({'SMILES': smiles, 'Tanimoto_Similarity': similarity_values}))
#
# # Save the similarity DataFrame to a CSV file
# similarity_df.to_csv('tanimoto_reactants_to_scaffold.csv', index=False)

# print('total number of reactions:', len(rxn_list))
# print('set of reactions:', len(set(rxn_list)))
# print(set(rxn_list))
# print(set(price_list))
df.to_csv('Libinvent_manifold_labelled.csv')
# Close the sqlitedict and the database connection
table.close()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        table_file_path = sys.argv[1]
        print(f'Cache path provided: {table_file_path}')
    else:
        print('No cache path provided, please provide path to cache file as an argument. '
              'Ex. python process_cache_output.py /path/to/cache/file')

    # Open the table as a sqlitedict
    table = sqlitedict.SqliteDict(
        filename=table_file_path,
        flag='c', autocommit=False, tablename='unnamed')

    # Only looking at 1 compound
    one_compound = False
    i = 0
    rxn_list = []
    df = pd.DataFrame()
    similarity_df = pd.DataFrame()
    for key, value in table.items():
        i += 1
        dictionary = value
        if i == 1:
            df = analyse_routes(dictionary)
        else:
            df = df.append(analyse_routes(dictionary))
        rxn_list, price_list = all_reactions_seen(dictionary, rxn_list)
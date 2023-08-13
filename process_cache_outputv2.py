import sqlitedict
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
import sys
import sqlite3
import os
import argparse
import csv


class PosteraTable:
    def __init__(self, table_file_path):
        self.table = sqlitedict.SqliteDict(filename=table_file_path, flag='c', autocommit=False, tablename='unnamed')
        self.sqlite3table = sqlite3.connect(table_file_path)

    def get_reactions(self):
        return self.table.items()

    def print_table(self):
        def truncate(value, length):
            return (value[:length] + '...') if len(value) > length else value

        # Set the length to truncate the value
        length_to_truncate = 100

        # Print all key-value pairs with truncated values
        for key, value in self.table.items():
            truncated_value = truncate(str(value), length_to_truncate)
            print("Key:", key)
            print("Value:", truncated_value)
            print()

    def close(self):
        self.table.close()

class SuperstructureAnalysis:
    def __init__(self, smiles_list: list(), results_dir, source_csv_name):
        self.query_smiles_list = smiles_list
        self.results_dir = results_dir
        self.reactions = []
        self.rxn_list = []
        self.price_list = []
        self.df = pd.DataFrame()
        self.source_csv_name = source_csv_name.split('.')[-2].split('/')[-1]

    def analyse_structures(self, query_smiles, dictionary):
        print(dictionary)
        return pd.DataFrame(data)

    def process_structures(self):
        for smile in self.query_smiles_list:
            value = superstructure_table.table.get(Chem.MolToSmiles(Chem.MolFromSmiles(smile)))
            if value is not None:
                dictionary = value
                self.df = pd.concat([self.df, self.analyse_structures(smile, dictionary)])
            else:
                print(f'SMILES string {self.query_smiles} not found in the database. Continuing to next SMILES string.')
                continue


class RouteAnalysis:
    def __init__(self, smiles_list: list(), results_dir, num_steps, source_csv_name):
        self.query_smiles_list = smiles_list
        self.results_dir = results_dir
        self.num_steps = num_steps
        self.reactions = []
        self.rxn_list = []
        self.price_list = []
        self.df = pd.DataFrame()
        self.source_csv_name = source_csv_name.split('.')[-2].split('/')[-1]

    def calculate_tanimoto_similarity(self, smiles1, smiles2):
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
        similarity = DataStructs.FingerprintSimilarity(fp1, fp2)
        print(similarity)
        return similarity

    def all_reactions_seen(self, dictionary):
        rxn_order = []
        price_list = []
        for route in dictionary['routes']:
            rxn_order = [reaction['name'] for reaction in route['reactions']]
            for molecule in route['molecules']:
                price_list.extend([entry.get('bbPriceRange') for entry in molecule['catalogEntries'] if molecule['isBuildingBlock']])
                price_list.extend([entry.get('scrPriceRange') for entry in molecule['catalogEntries'] if not molecule['isBuildingBlock']])

        self.rxn_list.extend(rxn_order)
        self.price_list.extend(price_list)

    def analyse_routes(self, smile: str(), dictionary):
        """
        Input: Single smiles string

        :param smile:
        :param dictionary:
        :return: DataFrame of organized metadata for each route found for query smiles.
        """
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
                # catalog_names.extend([entry.get('catalogName') for entry in route['molecules'] if entry['smiles'] in reactant_smiles])

            for molecule in route['molecules']:
                building_blocks.append(molecule['smiles'] if molecule['isBuildingBlock'] else None)
                catalog_names.extend(
                    [entry.get('catalogName') for entry in molecule['catalogEntries'] if molecule['isBuildingBlock']])
                lead_time.extend([entry.get('bbLeadTimeWeeks') for entry in molecule['catalogEntries'] if
                                  molecule['isBuildingBlock']])
            # building_blocks_for_rxn = [building_blocks[i] if i < len(building_blocks) else None for i in range(num_steps)]
            # catalog_name_for_building_block = [catalog_names[i] if i < len(catalog_names) else None for i in range(num_steps)]

            car_route = all(name in postera_terminology_reactions for name in encoded_reactions)
            # car_route = any(any(name in term for term in postera_terminology_reactions) for name in rxn_order)

            data.append({
                'SMILES': smile,
                #'CAR_route': car_route,
                'num_steps': num_steps,
                'rxn_order_first_to_last': rxn_order,
                'reactants': reactants,
                'BuildingBlocks': building_blocks,
                'catalogName_for_BuildingBlock': catalog_names,
            })

        return pd.DataFrame(data)

    def process_reactions(self):
        for smile in self.query_smiles_list:
            value = reaction_table.table.get(Chem.MolToSmiles(Chem.MolFromSmiles(smile)))
            if value is not None:
                dictionary = value
                self.df = pd.concat([self.df, self.analyse_routes(smile, dictionary)])
                self.all_reactions_seen(dictionary)
            else:
                print(f'SMILES string {self.query_smiles} not found in the database. Continuing to next SMILES string.')
                continue
        if self.num_steps is not None:
            # Saving csv of all reaction information to only include specified number of steps
            self.df = self.df[self.df['num_steps'] == self.num_steps]
            if self.num_steps == 1:
                print(
                    f'Saving all routes with {self.num_steps} step(s) to csv at {self.source_csv_name}_routes_{self.num_steps}_step.csv')
                self.df.to_csv(os.path.join(self.results_dir, f'{self.source_csv_name}_routes_{self.num_steps}_step.csv'), index=False)
            else:
                print(
                    f'Saving all routes with {self.num_steps} step(s) to csv at {self.source_csv_name}_routes_{self.num_steps}_steps.csv')
                self.df.to_csv(os.path.join(self.results_dir, f'{self.source_csv_name}_routes_{self.num_steps}_steps.csv'), index=False)
        else:
            print(f'Saving all routes for all compounds to csv at {self.source_csv_name}_routes.csv')
            self.df.to_csv(os.path.join(self.results_dir, f'{self.source_csv_name}_routes.csv'), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process cache output")
    parser.add_argument('-t',"--table_file_path", help="Path to the cache file", required=True)
    parser.add_argument('-r',"--results_dir", help="Directory for the results", required=True)
    parser.add_argument('-s',"--smiles", help="csv of SMILES you are looking for", required=True)
    parser.add_argument('-x', "--retrosynthesis", help="if performing a retrosynthesis search", action="store_true")
    parser.add_argument('-u', "--superstructure", help='if performing a superstructure search', action="store_true")
    parser.add_argument('-n',"--num_steps", type=int, help="Number of steps you are looking for")

    args = parser.parse_args()

    # TODO: Could parallelize search if searching for many SMILES
    # Load the smiles from the CSV file
    with open(args.smiles, 'r') as f:
        reader = csv.reader(f)
        next(reader, None)  # skip the headers
        smiles_list = list(reader)

    # Flatten the list if the CSV has only one column
    smiles_list = [item for sublist in smiles_list for item in sublist]

    print(f'Cache path provided: {args.table_file_path}')
    print(f'Results directory provided: {args.results_dir}')
    if args.retrosynthesis:
        print(f'---Performing retrosynthesis search---')
        print(f'Number of steps you are looking for: {args.num_steps}')
        reaction_table = PosteraTable(args.table_file_path)

        # Pass the list of smiles to the RouteAnalysis class
        route_analysis = RouteAnalysis(smiles_list=smiles_list, results_dir=args.results_dir, num_steps=args.num_steps, source_csv_name=args.smiles)
        route_analysis.process_reactions()
        reaction_table.close()

    elif args.superstructure:
        print(f'---Performing superstructure search---')
        superstructure_table = PosteraTable(args.table_file_path)

        # Pass the list of smiles to the SupertructureAnalysis class
        super_analysis = SuperstructureAnalysis(smiles_list=smiles_list, results_dir=args.results_dir, source_csv_name=args.smiles)
        super_analysis.process_structures()
        super_analysis.close()


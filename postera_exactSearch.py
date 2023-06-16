import numpy as np
from postera_base import Postera_base


class Postera_exactSearch(Postera_base):


    def __init__(self, verbose=True):
        
        
        url = "https://api.postera.ai/api/v1/exact/"
        batch_size = 100

        super(type(self), self).__init__(  url, batch_size, verbose=verbose, 
                                            cache_fname_tags="")
        
    @property
    def computation_type_tag(self):
        return "cost"
        
    def _retrieveLeadTimeAndPriceFromRecord(self, record):
        if isinstance(record, dict):
            record = record.get("catalogEntries",record)

        t= np.inf
        cost = np.inf
        for entry in record:
            if "leadTimeWeeks" in entry:
                t = min(t, entry["leadTimeWeeks"])
                up, lp= entry["upperPricePerMg"] , entry["lowerPricePerMg"]
                cost = min(cost, lp + 0.5*(up - lp) )
        return (t, cost)
                
    def geResultFromRecord(self, record):
        if 'error' in record:
            result = (np.nan, np.nan)
        else:
            result = self._retrieveLeadTimeAndPriceFromRecord(record)
        return result
        
    def geResultFromCache(self, smi):
        return self.geResultFromRecord(self.cache[smi])

    def _processResultForCache(self, record):
        return record


def test():
    record = {'results': [{'catalogName': 'emolecules', 'catalogId': '488000', 'smiles': 'CCCCOC', 'link': 'https://www.emolecules.com/search/#?querytype=emoleculesid&p=1&query=488000', 'inchikeyMatches': {'exact': True, 'parent': True, 'connectivity': True}, 'leadTimeWeeks': 4.0, 'lowerPricePerMg': 0.0, 'upperPricePerMg': 5.0}, {'catalogName': 'mcule', 'catalogId': 'MCULE-9814272675', 'smiles': 'O(C)CCCC', 'link': 'https://mcule.com/MCULE-9814272675', 'inchikeyMatches': {'exact': True, 'parent': True, 'connectivity': True}, 'leadTimeWeeks': 4.0, 'lowerPricePerMg': 0.0, 'upperPricePerMg': 5.0}, {'catalogName': 'mcule_ultimate', 'catalogId': 'CXBDYQVECUFKRK-UHFFFAOYSA-N', 'smiles': 'C(OC)CCC', 'link': 'https://ultimateapp.mcule.com/search/?t=exact&q=CXBDYQVECUFKRK-UHFFFAOYSA-N', 'inchikeyMatches': {'exact': True, 'parent': True, 'connectivity': True}, 'leadTimeWeeks': 5.0, 'lowerPricePerMg': 250.0001, 'upperPricePerMg': 10000.0}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL16584691', 'smiles': 'CCO.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL16584691/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL210725', 'smiles': 'CC(O)CO.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL210725/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL265948', 'smiles': 'N.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL265948/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL921985', 'smiles': 'CS(O)(=O)=O.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL921985/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL4943823', 'smiles': 'CO.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL4943823/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL5765992', 'smiles': 'CCCCOC.OCC(O)CO', 'link': 'https://www.surechembl.org/chemical/SCHEMBL5765992/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL9614568', 'smiles': '[Na].CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL9614568/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL9788284', 'smiles': 'CC(O)=O.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL9788284/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL10532183', 'smiles': 'COCCO.CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL10532183/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL10581975', 'smiles': 'CCCCOC.CCOC(C)=O', 'link': 'https://www.surechembl.org/chemical/SCHEMBL10581975/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'surechembl', 'catalogId': 'SCHEMBL11574005', 'smiles': '[Na].[Na].CCCCOC', 'link': 'https://www.surechembl.org/chemical/SCHEMBL11574005/', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '56846530', 'smiles': 'CCCC[OH+]C', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/56846530', 'inchikeyMatches': {'exact': True, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '139781945', 'smiles': 'CCCCOC.Cl.Cl', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/139781945', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '86649813', 'smiles': 'CCCCOC.CS(=O)(=O)O', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/86649813', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '12420856', 'smiles': '[CH2-]CCCOC.[CH2-]CCCOC.[Zn+2]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/12420856', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '149980020', 'smiles': '[2H]COCCCC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/149980020', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '87752651', 'smiles': 'CCCCOC.CO', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/87752651', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '23371081', 'smiles': 'CCCCOC.N', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/23371081', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '21523331', 'smiles': 'CCCCOC.[Br-]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/21523331', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '59535762', 'smiles': '[2H]C([2H])([2H])OCCCC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/59535762', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '11679839', 'smiles': '[CH2-]CCCOC.[Cl-].[Mg+2]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/11679839', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '10877055', 'smiles': '[CH2-]CCCOC.[Li+]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/10877055', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '88488712', 'smiles': 'CCCCOC.[Na]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/88488712', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '22026576', 'smiles': 'CCCCOC.OCC(O)CO', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/22026576', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '17844943', 'smiles': 'CCCCOC.CCOC(C)=O', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/17844943', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '57142162', 'smiles': '[3H]COCCCC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/57142162', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '117966804', 'smiles': 'CCCCOC.CCO', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/117966804', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '88516035', 'smiles': 'CC(=O)O.CCCCOC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/88516035', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '88625945', 'smiles': 'CCCCOC.COCCO', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/88625945', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '88800055', 'smiles': 'CCCCOC.[Na].[Na]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/88800055', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '18373627', 'smiles': 'CC(O)CO.CCCCOC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/18373627', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '15976975', 'smiles': '[CH2-]CCCOC.[Cl-].[Mg+2]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/15976975', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '134979088', 'smiles': '[Be+2].[CH2-]CCCOC.[CH2-]CCCOC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/134979088', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '12338', 'smiles': 'CCCCOC', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/12338', 'inchikeyMatches': {'exact': True, 'parent': True, 'connectivity': True}}, {'catalogName': 'pubchem', 'catalogId': '58249694', 'smiles': '[CH2-]CCCOC.[K+]', 'link': 'https://pubchem.ncbi.nlm.nih.gov/compound/58249694', 'inchikeyMatches': {'exact': False, 'parent': True, 'connectivity': True}}]}
    
    print( Postera_exactSearch( verbose=False)._retrieveLeadTimeAndPriceFromRecord(record.get("results", record)))
    
if __name__ == "__main__":

    '''

    '''
    
#    test(); import sys; sys.exit(0)
    mols = Postera_exactSearch.parser()
    psa= Postera_exactSearch( verbose=False)
    preds = psa.search_from_molecules_generator( mols )
    print(preds)
    
#    import sys, os
#    import pandas as pd
#    from math import isinf
#    df = pd.read_csv( os.path.expanduser(sys.argv[1]))
#    df["in_catalogues"] = list(map(lambda x: not isinf(x[0]), preds ))
#    df.to_csv( "/home/sanchezg/oxford/myProjects/spend_65K_jun/manifold/order/manifold_final_order_isInCatalogue.csv", index=False)
    
    
    '''
python postera_exactSearch.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
    '''
    

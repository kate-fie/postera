import requests
import numpy as np
import itertools
from postera_base import Postera_base


class Postera_similaritySearch(Postera_base):

    SAVE_STEPS_FREQ = 5

    def __init__(self, banned_catalogues=["pubchem"], vendors = None, patents=None, useThirdParties=False, withPurchaseInfo=False,
                       maxPrice=None, maxResultPages=5, verbose=True, *args, **kwargs):
        

        url = "https://api.postera.ai/api/v1/similarity/"
        batch_size = 1
        self.banned_catalogues = banned_catalogues
        self.vendors = vendors
        self.patents = patents
        self.useThirdParties = useThirdParties
        self.withPurchaseInfo = withPurchaseInfo
        self.maxPrice = maxPrice
        self.maxResultPages = maxResultPages

        super(type(self), self).__init__(  url, batch_size, verbose=verbose, 
                                            cache_fname_tags="", *args, **kwargs)
        
    @property
    def computation_type_tag(self):
        return "similarity"
                
    def geResultFromRecord(self, record):

        if 'error' in record:
            valid_found = None
        else:
            valid_found = []
            for result in record["results"]:
                validCatalogue = len([ cat for cat in result["catalogEntries"] if cat["catalogName"] not in self.banned_catalogues ]) >0
                if validCatalogue:
                    smi = result["smiles"]
                    valid_found.append( (result["similarity"], smi, result["catalogEntries"]))
                
        return valid_found
        
    def geResultFromCache(self, smi):
        return self.cache[smi]

    def _processResultForCache(self, record):
        return self.geResultFromRecord(record)

    def makeQuery(self, query_data, queryNum=-1, nextPage=1):
        if nextPage is None:
          return {}
        if nextPage >=1:
          query_data["page"] = nextPage
          
        r = requests.post(self.url, json=query_data, headers={"X-API-KEY": Postera_base.MANIFOLD_API_KEY})
        if r.ok:
          json_out = r.json()
          next_page = json_out["nextPage"]
          return json_out, next_page
        else:
          print( "query %d failed!!"%j )
          print( r )
          if r.status_code == 429:
              self.save_cache(closeDb=True)
              print("Too many queries!")
              os._exit(1)
                
    def search(self, mol_list):

        mol_list, searchResults,  smiles_list, n_to_compute = self.pre_search(mol_list)

        query_data = {}
        if self.vendors:
          query_data["vendors"] = self.vendors
        if self.patents:
          query_data["patentDatabases"] = self.patentDatabases
        else:
          query_data["patentDatabases"] = []
        if self.useThirdParties:
          query_data["queryThirdPartyServices"] = True
        if self.withPurchaseInfo:
            query_data["withPurchaseInfo"] = True
        if self.maxPrice:
          query_data["maxScrPricePerMg"] = self.maxPrice
          query_data["maxBbPricePerG"] = self.maxPrice
          
        for j, smi  in smiles_list:
          if self.verbose: print("lauching query for smiles number %d (%s)..."%(j, self.url))
          #vendors = None, patents=None, useThirdParties=
          query_data["smiles"] = smi
          
          list_jsonOut = []
          nextPage = 1
          while nextPage is not None:
            json_out, nextPage = self.makeQuery(query_data, queryNum=j, nextPage=nextPage)
            if self.verbose: print("Page: %d" % (nextPage))

            list_jsonOut.append( json_out)
            if self.maxResultPages and nextPage > self.maxResultPages:
                break
          json_out["results"] = list(itertools.chain.from_iterable( js["results"] for js in list_jsonOut if "results" in js))
          self.storeResultToCache(smi, json_out)
          searchResults[j] = self.geResultFromRecord(json_out)
  
          if (j+1) % type(self).SAVE_STEPS_FREQ == 0:
            self.save_cache()

        return list(zip(mol_list, searchResults))
        
        
if __name__ == "__main__":

    '''

    '''
    
    mols = Postera_similaritySearch.parser()
    # mols = Postera_similaritySearch.parser(["postera_similaritySearch.py", "~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv", "1", "10"])

#    psa= Postera_similaritySearch( verbose=False)
#     psa= Postera_similaritySearch(cache_fname="prueba_cache.sqlite", verbose=True, maxResultPages=5, maxPrice=100, withPurchaseInfo=True, vendors=["mcule", "mcule_ultimate",])
    preds = psa.search_from_molecules_generator( mols )
    print(preds)
    
#    import sys, os
#    import pandas as pd
#    from math import isinf
#    df = pd.read_csv( os.path.expanduser(sys.argv[1]))
#    df["in_catalogues"] = list(map(lambda x: not isinf(x[0]), preds ))
#    df.to_csv( "/home/sanchezg/oxford/myProjects/spend_65K_jun/manifold/order/manifold_final_order_isInCatalogue.csv", index=False)
    
    
    '''
python postera_similaritySearch.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
    '''
    

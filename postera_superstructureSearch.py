import requests
import numpy as np
import itertools
from postera_base import Postera_base
import os


class Postera_superstructureSearch(Postera_base):
    SAVE_STEPS_FREQ = 5

    def __init__(self, catalogues=["mcule_ultimate", "enamine_real", "wuxi_bb_screening", "sigma", "generic",
                                     "molport", "emolecules", "mcule", "wuxi_galaxi", "enamine_bb", "enamine_made"], cache_fname_tags=[], vendors=None, patents=None, useThirdParties=False,
                 withPurchaseInfo=True,
                 maxPrice=None, maxResultPages=3, verbose=True, *args, **kwargs):

        url = "https://api.postera.ai/api/v1/superstructure/"
        batch_size = 1
        self.catalogues = catalogues
        self.cache_fname_tags = cache_fname_tags
        self.vendors = vendors
        self.patents = patents
        self.useThirdParties = useThirdParties
        self.withPurchaseInfo = withPurchaseInfo
        self.maxPrice = maxPrice
        self.maxResultPages = maxResultPages
        self.patentDatabases = ["all"]

        super(type(self), self).__init__(url, batch_size, verbose=verbose,
                                         cache_fname_tags=cache_fname_tags, *args, **kwargs)

    @property
    def computation_type_tag(self):
        return "superstructure"

    def geResultFromRecord(self, record):

        if 'error' in record:
            valid_found = None
        else:
            valid_found = []
            for result in record["results"]:
                print(result["catalogEntries"])
                validCatalogue = len(
                    [cat for cat in result["catalogEntries"] if cat["catalogName"] in self.catalogues]) > 0
                if validCatalogue:
                    smi = result["smiles"]
                    valid_found.append((smi, result["catalogEntries"]))

        return valid_found

    def geResultFromCache(self, smi):
        return self.cache[smi]

    def _processResultForCache(self, record):
        return self.geResultFromRecord(record)

    def makeQuery(self, query_data, queryNum=-1, nextPage=1):
        """

        :param query_data:
        :param queryNum:
        :param nextPage: results are in pages, 1 is the first page of results, 2 is the second page of results, etc.
        :return:
        """
        if nextPage is None:
            return {}
        if nextPage >= 1:
            query_data["page"] = nextPage

        r = requests.post(self.url, json=query_data, headers={"X-API-KEY": Postera_base.MANIFOLD_API_KEY})
        if r.ok:
            json_out = r.json()
            next_page = json_out["nextPage"] # None when there is not another page of results
            return json_out, next_page
        else:
            print("query %d failed!!")
            print(r)
            if r.status_code == 429:
                self.save_cache(closeDb=True)
                print("Too many queries!")
                os._exit(1)

    def search(self, mol_list):

        mol_list, searchResults, smiles_list, n_to_compute = self.pre_search(mol_list) # Search in cache

        query_data = {}
        if self.vendors:
            query_data["vendors"] = self.vendors
        if self.patents:
            query_data["patentDatabases"] = self.patentDatabases
        else:
            query_data["patentDatabases"] = ["all"]
        if self.useThirdParties:
            query_data["queryThirdPartyServices"] = False
        if self.withPurchaseInfo:
            query_data["withPurchaseInfo"] = True
        if self.maxPrice:
            query_data["maxScrPricePerMg"] = self.maxPrice
            query_data["maxBbPricePerG"] = self.maxPrice

        for j, smi in smiles_list:
            if self.verbose: print("lauching query for smiles number %d (%s)..." % (j, self.url))
            # vendors = None, patents=None, useThirdParties=
            query_data["smiles"] = smi

            list_jsonOut = []
            nextPage = 1
            while nextPage is not None:
                if self.verbose: print("Page: %d" % (nextPage))
                if self.maxResultPages and nextPage > self.maxResultPages:
                    break
                json_out, nextPage = self.makeQuery(query_data, queryNum=j, nextPage=nextPage)

                list_jsonOut.append(json_out)
            json_out["results"] = list(
                itertools.chain.from_iterable(js["results"] for js in list_jsonOut if "results" in js))
            self.storeResultToCache(smi, json_out)
            searchResults[j] = self.geResultFromRecord(json_out)

            if (j + 1) % type(self).SAVE_STEPS_FREQ == 0:
                self.save_cache()

        return list(zip(mol_list, searchResults))


if __name__ == "__main__":
    '''

    '''

    mols = Postera_superstructureSearch.parser()
    psa = Postera_superstructureSearch(catalogues=["mcule_ultimate", "generic", "molport", "mcule", "enamine_bb"],
                                       cache_fname_tags='enamine')
    preds = psa.search_from_molecules_generator(mols)
    print(preds) # Want to print preds with query molecule and found molecules

    '''
python postera_superstructureSearch.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
    '''


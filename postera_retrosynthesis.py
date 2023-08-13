import numpy as np
import time
from postera_base import Postera_base


class Postera_RsPlan(Postera_base):

    def __init__(self, maxSearchDepth=4, use_only_diamond_reactions=False, ignore_zero_steps=False,
                 catalogues=["generic", "mcule"], customBuildingBlock=[], cache_fname_tags='', verbose=True):

        url = "https://api.postera.ai/api/v1/retrosynthesis/"
        batch_size = 10

        self.maxSearchDepth = maxSearchDepth
        self.ignore_zero_steps = ignore_zero_steps
        self.catalogues = catalogues

        cache_fname_tags = cache_fname_tags
        if use_only_diamond_reactions:
            cache_fname_tags += "_diamond"
        cache_fname_tags += "-" + "-".join(catalogues)

        self.use_only_diamond_reactions = use_only_diamond_reactions
        self.customBuildingBlock = customBuildingBlock

        super(type(self), self).__init__(url, batch_size, verbose=verbose,
                                         cache_fname_tags=cache_fname_tags)

    @property
    def computation_type_tag(self):
        return "retrosynthesis"

    def _getMinNSteps(self, result):

        nSteps = float("inf")
        if len(result) == 0 or "error" in result:
            return np.nan
        for route in result["routes"]:
            nSteps_ = len(route["reactions"])
            if self.ignore_zero_steps and nSteps_ == 0:
                continue
            nSteps = min(nSteps, nSteps_)
        return nSteps

    def geResultFromRecord(self, record):
        if 'error' in record:
            score = np.nan
        else:
            score = self._getMinNSteps(record)
        return score

    def geResultFromCache(self, smi):
        return self._getMinNSteps(self.cache[smi])

    def _processResultForCache(self, record):
        return record

    def prepateQueryData(self, selected_mols):
        time.sleep(2)
        url, data = super(type(self), self).prepateQueryData(selected_mols)

        data["maxSearchDepth"] = self.maxSearchDepth
        data["catalogs"] = self.catalogues
        if self.use_only_diamond_reactions:
            data["reactionTag"] = "diamond_robotic_synthesis"
        if len(self.customBuildingBlock) > 0:
            data["customBuildingBlocks"] = self.customBuildingBlock
        return url, data


if __name__ == "__main__":

    '''

    '''
    mols = Postera_RsPlan.parser()
    psa = Postera_RsPlan(maxSearchDepth=4, use_only_diamond_reactions=True,
                         ignore_zero_steps=True,
                         catalogues=["mcule_ultimate", "enamine_real", "wuxi_bb_screening", "sigma", "generic",
                                     "molport", "emolecules", "mcule", "wuxi_galaxi", "enamine_bb", "enamine_made"],
                         cache_fname_tags='scaffold_restrict',
                         verbose=True)
    preds = psa.search_from_molecules_generator(mols)
    print(preds)

    for smi, pred in zip(mols, preds):
        with open("Libinvent_retrosynthesis_output_June20.csv", "w") as f:
            f.write("%s\t%s\n" % (smi, pred))

    '''
python postera_retrosynthesis.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1 200
python postera_retrosynthesis.py  ['NC(=O)c1cccc(NC(=O)NCc2ccc(Cl)cc2)c1',
 'Cc1csc(CNC2CCNC2=O)n1',
 'Cc1cc2oc(=O)c3c(c2c2c1C(=O)CC(c1ccc4c(c1)OCO4)O2)CCCC3',
 'COc1cc(CNC(=O)NCc2ccc3c(c2)OCO3)ccc1O',
 'CC(=O)CCC(=O)C1CCC(C2=CCC3=C2CC2(C3)N=Nc3ccccc32)NC1']
    '''

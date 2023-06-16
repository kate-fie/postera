import numpy as np
from postera_base import Postera_base

class Postera_SA(Postera_base):


    def __init__(self, use_fast_instead=False, use_only_diamond_reactions=False,
                 return_number_steps_instead = False, verbose=True):
        #super(Postera_SA, self).__init__()
        cache_fname_tags=""
        if use_fast_instead:
            assert not return_number_steps_instead, "Error, ""minNumSteps"" only valid when use_fast_instead=False"
            url = "https://api.postera.ai/api/v1/synthetic-accessibility/fast-score/"
            self.retrieve_score_tag = "fastSAScore"
            batch_size = 100
            cache_fname_tags += "_fast"
        else:
            url = "https://api.postera.ai/api/v1/synthetic-accessibility/retrosynthesis/"
            self.retrieve_score_tag = "score"
            batch_size = 10 #20
            
        self.return_number_steps_instead = return_number_steps_instead
        self.minNumStepsTag = "minNumSteps"
        if self.return_number_steps_instead:
            assert not use_fast_instead, "Error, ""minNumSteps"" only valid when use_fast_instead=False"
            
        self.cache_idx = 1 if self.return_number_steps_instead else 0
        
        self.use_only_diamond_reactions = use_only_diamond_reactions
        if use_only_diamond_reactions:
            cache_fname_tags += "_diamond"

        super(type(self), self).__init__(  url, batch_size, verbose=verbose, 
                                            cache_fname_tags=cache_fname_tags)
        
    @property
    def computation_type_tag(self):
        return "SA"
        
    def _retrieveScoreAndStepsFromRecord(self, record):
        record = record.get("SAData",record)
        return (1 - record[self.retrieve_score_tag], 
                record.get(self.minNumStepsTag,None) )
                
    def geResultFromRecord(self, record):
        if 'error' in record:
            score = np.nan
        else:
            sa_score, stesps = self._retrieveScoreAndStepsFromRecord(record)
            if self.return_number_steps_instead:
                score = stesps
            else:
                score = sa_score
        return score
        
    def geResultFromCache(self, smi):
        return self.cache[smi][self.cache_idx]

    def _processResultForCache(self, record):
        if 'error' in record:
            sa_score = np.nan
            stesps = np.nan
        else:
            sa_score, stesps = self._retrieveScoreAndStepsFromRecord(record)
        return (sa_score, stesps )
        
    
    def prepateQueryData(self, selected_mols):
        url, data =  super(type(self),self).prepateQueryData(selected_mols)
        if self.use_only_diamond_reactions:
             data["reactionTag"]= "diamond_robotic_synthesis"
        return url, data
        
       
if __name__ == "__main__":

    '''

    '''
    mols = Postera_SA.parser()
    psa= Postera_SA(use_fast_instead=False, use_only_diamond_reactions=False,
                 return_number_steps_instead=False, verbose=True)
    preds = psa.search_from_molecules_generator(mols)
    print(preds)

#    print(preds)
    
####    import pandas as pd
####    df = pd.read_csv("/home/sanchezg/oxford/myProjects/spend_65K_jun/manifold/order/manifold_final_order_isInCatalogue.csv")
####    df["sa_robot"] = preds
####    df.to_csv("manifold_final_order_isInCatalogue_robotSA.csv", index=False)
    

    '''
python postera_sa.py  ~/oxford/myProjects/deepLearningCompounds/rawData/compoundsCost/Mcule/_mcule_400K_testSet.csv 1
    '''
    

import sys, os, shutil
import requests
import gzip
from subprocess import check_output
import atexit
import tempfile
import tqdm
import numpy as np

from itertools import chain, islice
from rdkit import Chem
from sqlitedict import SqliteDict 
from abc import ABC, abstractmethod



class Postera_base(ABC):

    MANIFOLD_API_KEY = os.environ["MANIFOLD_API_KEY"]
    assert MANIFOLD_API_KEY, "Error, you need to set the environmental MANIFOLD_API_KEY"
    SAVE_STEPS_FREQ = 50

    def __init__(self, url, batch_size, verbose=True, cache_fname_tags="", 
                    cache_fname= os.path.abspath(os.path.join(__file__, "../cache/cache_postera_%(computation_type_tag)s%(cache_fname_tags)s.sqlite"))):
        

        self.url = url
        self.batch_size = batch_size
        computation_type_tag = self.computation_type_tag
        self.cache_fname = cache_fname%locals()
        self.verbose = verbose
        self.cache_dirty=False
        
        self.tmpdir = tempfile.TemporaryDirectory()
        self.cache_working_name = os.path.join(self.tmpdir.name, os.path.basename(self.cache_fname))
        if os.path.isfile(self.cache_fname):
            shutil.copyfile(self.cache_fname, self.cache_working_name)
        self.cache = SqliteDict(self.cache_working_name, autocommit=False)
        if self.verbose: print("Number of entries in cache:", len(self.cache))
        atexit.register( self.save_cache, closeDb=True )

    
    @property
    @abstractmethod
    def computation_type_tag(self):
        raise NotImplementedError() 
        
    @abstractmethod
    def geResultFromRecord(self, record):
        raise NotImplementedError()
        
        
    @abstractmethod
    def geResultFromCache(self, smi):
#        return self.cache[smi]
        raise NotImplementedError()

    @abstractmethod
    def _processResultForCache(self, record):
        raise NotImplementedError()
        
    def storeResultToCache(self, smi, record):
        self.cache[smi] = self._processResultForCache(record)
        self.cache_dirty=True
    
    def prepateQueryData(self, selected_mols_smis):
        data = {}
        if len(selected_mols_smis) == 1:
            url = self.url
            data["smiles"] = selected_mols_smis[0][1]
        else:
            url = self.url+"batch/"
            data["smilesList"] = [ x[1] for x in selected_mols_smis]
        return url, data
        
        
    def pre_search(self, mol_list):
        mol_list = list(mol_list)
        n_mols = len(mol_list)
        scores = [None]*n_mols
        smiles_list = []

        for i, mol in enumerate(mol_list):
            if mol is not None:
                smi = Chem.MolToSmiles(mol)
                if smi in self.cache:
                    scores[i] = self.geResultFromCache(smi)
                else:
                    smiles_list.append( (i, smi) )
            else:
                scores[i] = None

        n_to_compute = len(smiles_list)
        return mol_list, scores,  smiles_list, n_to_compute
        
    def search(self, mol_list):

        mol_list, scores,  smiles_list, n_to_compute = self.pre_search(mol_list)

        i = 0
        for i in range( n_to_compute//self.batch_size +int(bool(n_to_compute%self.batch_size))):
          selected_mols_smis = smiles_list[i*self.batch_size:(i+1)*self.batch_size]
          url, data = self.prepateQueryData( selected_mols_smis)
          if self.verbose: print("lauching query %d (%s)..."%(i, url))
          r = requests.post(url, json=data, headers={"X-API-KEY": Postera_base.MANIFOLD_API_KEY})
          if r.ok:
            json_out = r.json()
            if len(selected_mols_smis) >1 :
                try:
                    results = json_out["results"]
                except KeyError:
                    print(data)
                    print(url)
                    print( json_out)
                    raise
            else:
                results = [ json_out.get("results", json_out)]
            for (j, smi), result in zip(selected_mols_smis, results):
                return_value = self.geResultFromRecord(result)
                self.storeResultToCache(smi, result)
                scores[j] = return_value

          else:
            print( "query %d failed!!"%i )
            print( r )
            if r.status_code == 429:
                self.save_cache(closeDb=True)
#                raise Exception("Too many queries")
                print("Too many queries!")
                sys.exit(1)
                
                
            for (j, smi) in selected_mols_smis:
                scores[j] = np.nan
          if (i+1) % Postera_base.SAVE_STEPS_FREQ == 0:
            self.save_cache()

        return list(zip(mol_list, scores))
            
    def save_cache(self, closeDb=False):
        if self.cache_dirty:
            if self.verbose: print("saving cache", self.cache_working_name, self.cache_fname)
            self.cache.commit()
            if closeDb:
                self.cache.close()
                os.rename(self.cache_working_name, self.cache_fname)
                del self.tmpdir
            else:
                shutil.copyfile(self.cache_working_name, self.cache_fname)
        elif closeDb:
            os.rename(self.cache_working_name, self.cache_fname)
            del self.tmpdir

    @classmethod
    def _batcher(cls, iterable, batch_size):
        return iter(lambda: tuple(islice(iter(iterable), batch_size)), ())
    
    def search_from_molecules_generator(self, mol_gen ):
                          
        super_batch_size = self.batch_size * 10 # self.batch_size*4
        
        mol_results = chain.from_iterable( ( self.search( list(batch) ) for batch in self._batcher(mol_gen, super_batch_size ) ) ) 
        
        smiles_results = map(lambda m_r: ( Chem.MolToSmiles(m_r[0]), m_r[1]), filter(lambda m_r: m_r[0] is not None, mol_results ))
        try:
            smiles, preds = zip(* smiles_results )
        except ValueError:
            return None
        return preds
        

    @classmethod
    def parser(cls, argv= sys.argv):
        fname_or_smis = os.path.expanduser(argv[1])
        if os.path.isfile(fname_or_smis):
            if fname_or_smis.endswith(".sdf"):
                suppl = Chem.SDMolSupplier(fname_or_smis)
            else:
                smilesColumn=1
                if len(argv)>2:
                    smilesColumn = int(argv[2])
                else:
                    nameColumn=False
                    
                if fname_or_smis.endswith(".gz"):

                    class MySupplier():
                        def __init__(self, fname, delimiter=","):
                            self.fname = fname
                            self.delimiter = ","
                            self._len = int(check_output(["wc", "-l",fname_or_smis]).decode("utf-8").split()[0]) - 1
                        def __iter__(self):
                            with gzip.open(self.fname) as f:
                                header = f.readline()
                                for line in f:
                                    line = line.decode("utf-8")
#                                    print(line); print(repr(line.strip().split(self.delimiter)[smilesColumn]))#; input("enter")
                                    yield Chem.MolFromSmiles(line.strip().split(self.delimiter)[smilesColumn])
                        def __len__(self):
                            return self._len
                            
                    suppl = MySupplier(fname_or_smis)
                else:
                    if fname_or_smis.endswith(".csv"):
                        delimiter=","
                        titleLine=True
                    else:
                        delimiter="\t"
                        titleLine=False
                    suppl = Chem.SmilesMolSupplier(fname_or_smis, smilesColumn=smilesColumn, nameColumn=False, titleLine=titleLine, delimiter=delimiter)
            mols = ( mol for mol in tqdm.tqdm(suppl, total=len(suppl)) )
        else:
            smiles = fname_or_smis.split(",")
            mols = ( Chem.MolFromSmiles(smi) for smi in tqdm.tqdm(smiles, total=len(smiles)) )
            
        if len(argv)>4:
            take_first_n = int(argv[3])
            mols = islice(mols, take_first_n)
        return mols
   
   
        

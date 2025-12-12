###
#
# Code to process the data in annotationDB and process it to the
# sample data. Currently pulls from AnnotationGx
#
# Author: James Bannon
# Github: jbannon
#####


from typing import List, Dict
import pandas as pd
import numpy as np
import requests
import json
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import sys
from damply import dirs
import tqdm 
from itertools import product


bioassay_res = defaultdict(list)
for cpd in cids:
	assay_subset = all_bioassays[cpd]
	cpd_results = num_assays*['Not Measured']
	for assay in assay_subset:
		assay_idx = aid_to_idx[assay['aid']]
		outcome = 'Active' if assay['activity_outcome_method']==2 else 'Inactive'
		cpd_results[assay_idx]=outcome
	bioassay_res[cpd]=cpd_results


bioassay_res = pd.DataFrame(bioassay_res,index=[f"AID_{aid}" for aid in seen_bioassays])
bioassay_res.reset_index(drop=False, inplace= True, names = 'Assay')
bioassay_res.to_csv(dirs.PROCDATA / "sample_data" / "bioassays.csv",index = False)

targets = pd.unique(bdb['Target Name'])
n_targets = len(targets)
target_to_idx = {targets[i]:i for i in range(n_targets)}
binding_res = defaultdict(list)


for cpd in cids:
	cpd_results = n_targets*['Not Measured']
	cpd_binding = bdb[bdb['PubChem CID']==cpd]
	cpd_tgts = cpd_binding['Target Name'].values
	for tgt in cpd_tgts:
		_tgt_binding = cpd_binding[cpd_binding['Target Name']==tgt]
		if _tgt_binding['Ki (nM)'].values[0] == '':
			affinity = _tgt_binding['Kd (nM)'].values[0]
		else:
			affinity = _tgt_binding['Ki (nM)'].values[0]

		cpd_results[target_to_idx[tgt]]=affinity

	binding_res[cpd] = cpd_results

# print(binding_res)



bdb_res = pd.DataFrame(binding_res,index = targets)
bdb_res.reset_index(drop=False,inplace=True,names = 'Target')
bdb_res.to_csv(dirs.PROCDATA / "sample_data"/"binding_db.csv",index=False)
sys.exit()



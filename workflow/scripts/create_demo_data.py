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


def make_fingerprint_generators(
	radius_list:List[int]=[2,3],
	size_list:List[int]=[1024,2048]
	)->Dict[str,rdkit.Chem.rdFingerprintGenerator]:

	generators = {f"Morgan({r},{s})":AllChem.GetMorganGenerator(radius=r, fpSize = s) for r,s in product(radius_list,size_list)}
	return generators



bdb = pd.read_csv(dirs.PROCDATA / "bdb_temp.csv")
bdb['PubChem CID'] = [int(x) for x in bdb['PubChem CID'].values]
drug_universe = "https://annotationdb.bhklab.ca/compound/all"
jump_cp_compounds = pd.read_csv(dirs.RAWDATA / "JUMP-CP"/ "jump_cp_compounds.csv")
lincs_compounds = pd.read_csv(dirs.RAWDATA / "LINCS"/ "compoundinfo_beta.txt",sep="\t")
bbbp = pd.read_csv(dirs.RAWDATA / "BBBP.csv")


response = requests.get(drug_universe).json()
colData, all_bioassays = defaultdict(list),defaultdict(list)
 
fingerprint_generators = make_fingerprint_generators()
fp_data = {k:defaultdict(list) for k in fingerprint_generators.keys()}
errs = []
seen_bioassays, cids = [],[]
num_cpds_seen = 0


for drug_info in tqdm.tqdm(response):

	try:
		drug_query_url = f'https://annotationdb.bhklab.ca/compound/many?compounds={drug_info['cid']}&format=json'
		drug_details = requests.get(drug_query_url).json()[0]
	

		# Basic Drug Info
		colData['Compoound Name'].append(drug_info['name'])
		colData['PubChem CID'].append(drug_info['cid'])
		
		colData['SMILES'].append(drug_info['smiles'])
		colData['InChIKey'].append(drug_info['inchikey'])
		smiles_str = drug_info['smile']
		cid = drug_info['cid']
		cids.append(cid)
		

		# Query AnnotationDB for Drug Specifics
		
		
		colData['Molecular Formula'].append(drug_details['molecular_formula'])
		
		# Additional names/ids/descriptions
		colData['IUPAC Name'].append(drug_details['iupac_name'])
		colData['ChEMBL ID'].append(drug_details['molecule_chembl_id'])
		colData['PubChem 2D Fingerprint'].append(drug_details['fingerprint_2d'])
		

		# MOA and approval
		mechanisms = drug_details['mechanisms']
		if len(mechanisms)==0:
			colData['Mechanism of Action'].append("None")
		else:
			colData['Mechanism of Action'].append(drug_details['mechanisms'][0]['mechanism_of_action'])
		colData['FDA Approved'].append(drug_details['fda_approval'])
		

		## Add in Molecular Information (molecular weight + Lipinski Filters)
		## 	for the curious: https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five
		colData['Molecular Weight'].append(drug_details['molecular_weight'])
		colData['XlogP'].append(drug_details['xlogp'])
		colData['Hydrogen Bond Donors'].append(drug_details['h_bond_donor_count'])
		colData['Hydrogen Bond Acceptors'].append(drug_details['h_bond_acceptor_count'])
		colData['Exact Molecular Mass'].append(drug_details['exact_mass'])
		colData['DILI Severity'].append(drug_details['toxicity'][ "dili_severity_grade"])
		colData['DILI Annotation'].append(drug_details['toxicity']["dili_annotation"])
		colData['Hepatotoxiciy Likelihood (Detailed)'].append(drug_details['toxicity']["hepatotoxicity_likelihood_score"])
		score = drug_details['toxicity']["hepatotoxicity_likelihood_score"].split(":")[1].lstrip().split()[0]
		colData['Hepatotoxiciy Likelihood (Score)'].append(score)

		## Check Against The Broad Data
		l1k_subset = lincs_compounds[lincs_compounds['inchi_key']==drug_info['inchikey']]
		jump_subset = jump_cp_compounds[jump_cp_compounds['Metadata_InChIKey']==drug_info['inchikey']]
		
		if l1k_subset.shape[0]==0:
			colData['In L1000'].append(False)
			colData['L1000 ID'].append("-")
		else:
			colData['In L1000'].append(True)
			colData['L1000 ID'].append(l1k_subset['pert_id'].values[0])
		

		if jump_subset.shape[0]==0:
			colData['In JUMP-CP'].append(False)
			colData['JUMP-CP ID'].append("-")
		else:
			colData['In JUMP-CP'].append(True)
			colData['JUMP-CP ID'].append(jump_subset['Metadata_JCP2022'].values[0])
		

		blood_brain = bbbp[bbbp['cid']==cid]
		if blood_brain.shape[0]==0:
			colData['BBB Permeable'].append('Unknown')
		else:
			colData['BBB Permeable'].append(blood_brain['p_np'].values[0])

		# Make the fingerpints
		rdk_mol = Chem.MolFromSmiles(smiles_str)
		for fp_type in fingerprint_generators.keys():
			generator = fingerprint_generators[fp_type]
			fp = generator.GetCountFingerprintAsNumPy(rdk_mol)
			fp_data[fp_type]['CID'].append(cid)
			for dim in range(fp.shape[0]):
				fp_data[fp_type][f'V{dim+1}'].append(fp[dim].item())
			
		
	

		# Cache bioassays for post-processing
		all_bioassays[drug_info['cid']] = drug_details['bioassays']
		seen_bioassays.extend([assay['aid'] for assay in drug_details['bioassays']])
	except:
		errs.append(drug_info['cid'])



colData = pd.DataFrame(colData)

colData.to_csv(dirs.PROCDATA / "sample_data" / "colData.csv",index=False)
for fp_type in fp_data:
	fp_df = pd.DataFrame(fp_data[fp_type])
	fp_df.to_csv(dirs.PROCDATA / "sample_data" / f"{fp_type}.csv",index = False)



seen_bioassays = sorted(list(set(seen_bioassays)))
aid_to_idx = {seen_bioassays[i]:i for i in range(len(seen_bioassays))}
num_assays = len(seen_bioassays)
num_cpds = len(cids)
cid_to_idx = {cids[i]:i for i in range(len(cids))}


# bioassay_res = np.empty((num_assays,num_cpds))*np.nan


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



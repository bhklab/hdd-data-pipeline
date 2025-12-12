from typing import List, Dict, Union
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
from utils import process_single_drug
import argparse 






def make_fingerprint_generators(
	
	radius_list: List[int],

	dim_list: List[int]
	
	) -> Dict[str,rdkit.Chem.rdFingerprintGenerator]:

	generators = {f"Morgan({r},{s})":AllChem.GetMorganGenerator(radius=r, fpSize = s) for r,s in product(radius_list,size_list)}
	return generators




def main(
	db_url:str, 

	radius_list:List[int],

	dim_list:List[int],
	
	lincs_file:str,	

	jump_cp_file:str,	
	) -> None:

	# set up storage for results
	# default dict allows colData to be easily converted into a dataframe
	
	colData, all_bioassays = defaultdict(list),defaultdict(list)
	seen_bioassays, cids = [],[]
	
	error_cids = []
	
	# download all the compounds in annotationdb
	response = requests.get(db_url).json()
	
	# make the fingerprint generators
	fingerprint_generators = make_fingerprint_generators(radius_list, dim_list)
	
	# store as a dictionary for breaking out as files later
	fp_data = {k:defaultdict(list) for k in fingerprint_generators.keys()}


	jump_cp_compounds = pd.read_csv(jump_cp_file)
	lincs_compounds = pd.read_csv(lincs_file,sep="\t")

	
	for drug_info in response:
		

		try:
			drug_query_url = f'https://annotationdb.bhklab.ca/compound/many?compounds={drug_info['cid']}&format=json'
			drug_details = requests.get(drug_query_url).json()[0]

			utils.process_single_drug(
				drug_details = drug_details,
				colData = colData,
				fp_data = fp_data,
				all_bioassays = all_bioassays,
				seen_bioassays = seen_bioassays,
				lincs_compounds = lincs_compounds,
				jump_cp_compounds = jump_cp_compounds,
				fingerprint_generators = fingerprint_generators)
		
		except:
			errs.append(drug_info['cid'])
	

	# write and store colData
	colData = pd.DataFrame(colData)
	colData.to_csv(dirs.PROCDATA / "colData.csv",index=False)
	

	# process the bioassays since we have the data here.
	seen_bioassays = sorted(list(set(seen_bioassays)))
	aid_to_idx = {seen_bioassays[i]:i for i in range(len(seen_bioassays))}
	num_assays = len(seen_bioassays)
	num_cpds = len(cids)
	cid_to_idx = {cids[i]:i for i in range(len(cids))}

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
	bioassay_res.to_csv(dirs.PROCDATA /  "bioassays.csv",index = False)

	




if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog = 'make_ColData ',
		description = "Generate the colData and parse other annotationdb data")
	parser.add_argument('-u', help = "route to all compounds in database")
	parser.add_argument('-r', help = "radius list")
	parser.add_argument('-d', help = "dimension list")
	parser.add_argument('-l', help = "lincs")
	parser.add_argument('-j',help = 'jump cp')
	

	args = parser.parse_args()

	main(
		db_url = args.u,
		radius_list = args.r,
		dim_list = args.d,
		lincs_file = args.l,
		jump_cp_file = args.j)

	

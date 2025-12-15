import argparse
from collections import defaultdict
from itertools import product
from typing import Dict, List

import pandas as pd
import requests
import utils
from damply import dirs
from rdkit.Chem import AllChem
import tqdm


def main(
	db_url:str, 

	lincs_file:str,	

	jump_cp_file:str,	

	bbbp_file:str
	) -> None:

	# set up storage for results
	# default dict allows colData to be easily converted into a dataframe
	
	colData, all_bioassays = defaultdict(list),defaultdict(list)
	seen_bioassays, cids = [],[]
	
	error_cids = []
	
	# download all the compounds in annotationdb
	response = requests.get(db_url).json()


	# read in jumpcp and lincs
	jump_cp_compounds = pd.read_csv(jump_cp_file)
	lincs_compounds = pd.read_csv(lincs_file,sep="\t")
	blood_brain_perm = pd.read_csv(bbbp_file)

	
	for drug_info in tqdm.tqdm(response[:100]):
		

		try:

			drug_query_url = f'https://annotationdb.bhklab.ca/compound/many?compounds={drug_info['cid']}&format=json'
			drug_details = requests.get(drug_query_url).json()[0]

			utils.process_single_drug(
				drug_info,
				drug_details = drug_details,
				colData = colData,
				all_bioassays = all_bioassays,
				seen_bioassays = seen_bioassays,
				lincs_compounds = lincs_compounds,
				jump_cp_compounds = jump_cp_compounds,
				blood_brain_perm = blood_brain_perm
				)
		
		except:
			error_cids.append(drug_info['cid'])
	

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


	outpath = Path(dirs.PROCDATA / "experiments")
	bioassay_res = pd.DataFrame(bioassay_res,index=[f"AID_{aid}" for aid in seen_bioassays])
	bioassay_res.reset_index(drop=False, inplace= True, names = 'Assay')
	bioassay_res.to_csv(dirs.PROCDATA /  "bioassays.csv",index = False)

	




if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog = 'make_ColData ',
		description = "Generate the colData and parse other annotationdb data")
	parser.add_argument('-u', help = "route to all compounds in database")
	parser.add_argument('-l', help = "lincs")
	parser.add_argument('-j',help = 'jump cp')
	parser.add_argument('-b',help = 'blood brain barrier')
	

	args = parser.parse_args()

	main(
		db_url = args.u,
		lincs_file = args.l,
		jump_cp_file = args.j,
		bbbp_file = args.b)

	

from collections import defaultdict
from damply import dirs 
from rdkit import Chem
from typing import List, Dict
from itertools import product
from rdkit.Chem import AllChem
import pandas as pd  
from pathlib import Path

def make_fingerprint_generators(
	radius_list:List[int],
	dimension_list: List[int]
	) -> Dict[str, rdkit.Chem.rdFingerprintGenerator]:
	
	fingerprint_generators = {}
	for radius,dimension in product(radius_list,dimension_list):
		fpgen =  AllChem.GetMorganGenerator(radius=radius,fpSize = dimension)
		fingerprint_generators[f"Morgan({radius},{dimension})"] = fpgen

	return fingerprint_generators






def main(
	radius_list = [2,3],
	dimension_list = [512,1024,2048]
	):
	
	fingerprint_generators = make_fingerprint_generators(radius_list,dimension_list)
	fp_data = {k:defaultdict(list) for k in fingerprint_generators}
	# fp_data = defaultdict(list)
	colData = pd.read_csv(dirs.PROCDATA / "colData.csv",usecols = ['Pubchem CID', 'SMILES'])
	
	outpath = dirs.PROCDATA / "experiments"/ "fingerprints"
	Path(outpath).mkdir(parents=True,exist_ok=True)
	
	# Make the fingerpints
	for idx, row in colData.iterrows():
		cid = row['Pubchem CID']
		smiles_str = row['SMILES']
		
		rdk_mol = Chem.MolFromSmiles(smiles_str)
		

		for fp_type in fingerprint_generators:
			fp_data[fp_type]['CID'].append(cid)
			generator = fingerprint_generators[fp_type]
			fp = generator.GetCountFingerprintAsNumPy(rdk_mol)
			for j in  range(fp.shape[0]):
				vname = f"V{j+1}"
				fp_data[fp_type][vname].append(fp[j])
	
	for fp_type in fp_data:
		fp_str = ".".join(fp_type.replace("(",".").replace(")","").split(","))
		
		fp_matrix = pd.DataFrame(fp_data[fp_type])
		fp_matrix = fp_matrix.transpose()
		
		fp_matrix = fp_matrix.rename(columns = fp_matrix.iloc[0])

		fp_matrix = fp_matrix.iloc[1:,]
		fp_matrix = fp_matrix.reset_index(drop=True)
		fp_matrix.to_csv(outpath /f"{fp_str}.csv",index=False)
		

if __name__ == '__main__':
	# parser.add_argument('-r', help = "radius list")
	# parser.add_argument('-d', help = "dimension list")main()
	main()

	# radius_list = args.r,
	# dim_list = args.d,	
	# make the fingerprint generators
	# fingerprint_generators = make_fingerprint_generators(radius_list, dim_list)
	
	# store as a dictionary for breaking out as files later
	

	
		
	
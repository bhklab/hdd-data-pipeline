from collections import defaultdict
from damply import dirs 
from rdkit import Chem
from typing import List, Dict
from itertools import product
from rdkit.Chem import AllChem
import pandas as pd  


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
	fp_data = defaultdict(list)
	colData = pd.read_csv(dirs.PROCDATA / "colData.csv",usecols = ['Pubchem CID', 'SMILES'])

	for idx, row in colData.iterrows():
		cid = row['Pubchem CID']
		smiles_str = row['SMILES']
	# Make the fingerpints
		rdk_mol = Chem.MolFromSmiles(smiles_str)
		fp_data['CID'].append(cid)

		for fp_type in fingerprint_generators:
			generator = fingerprint_generators[fp_type]
			fp = generator.GetCountFingerprintAsNumPy(rdk_mol)
			fp_data[fp_type].append(fp)


	fp_data = pd.DataFrame(fp_data)
	fp_data.to_csv(dirs.PROCDATA / "experiments" / "fingerprints.csv")
		

if __name__ == '__main__':
	# parser.add_argument('-r', help = "radius list")
	# parser.add_argument('-d', help = "dimension list")main()
	main()

	# radius_list = args.r,
	# dim_list = args.d,	
	# make the fingerprint generators
	# fingerprint_generators = make_fingerprint_generators(radius_list, dim_list)
	
	# store as a dictionary for breaking out as files later
	

	
		
	
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from damply import dirs


def is_numeric(x):
	numeric_type = isinstance(x, int) or isinstance(x, float)
	return numeric_type


def main(coldata_path: str, bindingdb_path: str) -> None:
	colData = pd.read_csv(coldata_path, usecols=['Pubchem CID', 'SMILES'])
	binding_db = pd.read_csv(bindingdb_path)
	binding_db = binding_db.dropna(subset=['Ki (nM)', 'Kd (nM)'], how='all')
	binding_db = binding_db[binding_db['PubChem CID'].isin(colData['Pubchem CID'])]
	targets = list(pd.unique(binding_db['Target Name']))
	cpds = list(pd.unique(binding_db['PubChem CID']))

	num_cpds = len(cpds)
	num_tgts = len(targets)

	target_to_idx = {targets[i]: i for i in range(num_tgts)}
	seen_targets = []
	res = defaultdict(list)
	for cpd in cpds:
		cpd_data = num_tgts * [np.nan]
		# print(len(cpd_data))

		temp = binding_db[binding_db['PubChem CID'] == cpd]
		for idx, row in temp.iterrows():
			seen_targets.append(row['Target Name'])
			aff_col = 'Ki (nM)' if is_numeric(row['Kd (nM)']) else 'Kd (nM)'
			aff = row[aff_col]
			if isinstance(aff, str):
				if '>' in aff or '<' in aff:
					continue
				else:
					aff = float(aff)
			cpd_data[target_to_idx[row['Target Name']]] = aff
		res[cpd] = cpd_data
	# print(set(seen_targets))
	res = pd.DataFrame(res, index=targets)
	res = res.dropna(axis=1, how='all')
	res = res.dropna(axis=0, how='all')
	res.to_csv(dirs.PROCDATA / 'experiments' / 'binding_db.csv')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog='make_bdb_experiments',
		description='Generate BindingDB experiments matrix',
	)
	parser.add_argument(
		'-c',
		'--coldata',
		default=str(dirs.PROCDATA / 'colData.csv'),
		help='colData CSV path',
	)
	parser.add_argument(
		'-b',
		'--bindingdb',
		default=str(dirs.RAWDATA / 'BINDING_DB' / 'BindingDB_All_202512_cleaned.csv'),
		help='BindingDB cleaned CSV path',
	)
	args = parser.parse_args()

	main(coldata_path=args.coldata, bindingdb_path=args.bindingdb)

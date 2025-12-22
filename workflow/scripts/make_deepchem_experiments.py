import argparse
import pandas as pd
from damply import dirs


def process_deepchem_data(
	colData: pd.DataFrame,
	data: pd.DataFrame,
	index_name: str,
	convert_to_int: bool = False,
) -> pd.DataFrame:
	measurement_cols = [
		col_name for col_name in data.columns if col_name not in ['smiles', 'mol_id']
	]
	if convert_to_int:
		data[measurement_cols] = data[measurement_cols].astype(int)

	data = pd.merge(colData, data, left_on='SMILES', right_on='smiles')

	data = data[['Pubchem CID'] + measurement_cols].transpose()
	data = data.rename(columns=data.iloc[0])
	data = data.iloc[1:,]
	data = data.reset_index(names=index_name, drop=False)
	return data


def main(coldata_path: str, deepchem_subdir: str) -> None:
	colData = pd.read_csv(coldata_path, usecols=['SMILES', 'Pubchem CID'])
	deepchem_root = dirs.PROCDATA / deepchem_subdir
	clintox = pd.read_csv(deepchem_root / 'clintox.csv')
	tox21 = pd.read_csv(deepchem_root / 'tox21.csv')
	toxcast = pd.read_csv(deepchem_root / 'toxcast.csv')
	sider = pd.read_csv(deepchem_root / 'sider.csv')

	clintox = process_deepchem_data(colData, clintox, 'Clinical Tox Result')
	sider = process_deepchem_data(colData, sider, 'Side Effect')
	tox21 = process_deepchem_data(colData, tox21, 'Tox Assay', True)
	toxcast = process_deepchem_data(colData, toxcast, 'Tox Assay')

	clintox.to_csv(dirs.PROCDATA / 'experiments' / 'clintox.csv', index=False)
	sider.to_csv(dirs.PROCDATA / 'experiments' / 'sider.csv', index=False)
	tox21.to_csv(dirs.PROCDATA / 'experiments' / 'tox21.csv', index=False)
	toxcast.to_csv(dirs.PROCDATA / 'experiments' / 'toxcast.csv', index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog='make_deepchem_experiments',
		description='Generate DeepChem experiments matrix',
	)
	parser.add_argument(
		'-c',
		'--coldata',
		default=str(dirs.PROCDATA / 'colData.csv'),
		help='colData CSV path',
	)
	parser.add_argument(
		'-s',
		'--subdir',
		default='DEEP_CHEM',
		help='DeepChem subdirectory under procdata',
	)
	args = parser.parse_args()

	main(coldata_path=args.coldata, deepchem_subdir=args.subdir)

import pandas as pd 
from damply import dirs
import argparse 

def process_deepchem_data(
	colData: pd.DataFrame, 
	data: pd.DataFrame,
	index_name:str 
	) -> pd.DataFrame:
	
	measurement_cols = [col_name for col_name in data.columns if col_name not in ['smiles','mol_id']]
	data = pd.merge(colData,data,left_on = 'SMILES', right_on = 'smiles')

	data = data[['Pubchem CID'] + measurement_cols].transpose()
	data = data.rename(columns = data.iloc[0])
	data = data.iloc[1:,]
	data = data.reset_index(names = index_name,drop=False)
	return data
	 
def main():
	colData = pd.read_csv(dirs.PROCDATA / "colData.csv", usecols = ['SMILES','Pubchem CID'])
	clintox = pd.read_csv(dirs.PROCDATA / "DEEP_CHEM" /"clintox.csv")
	tox21 = pd.read_csv(dirs.PROCDATA / "DEEP_CHEM" / "tox21.csv")
	toxcast = pd.read_csv(dirs.PROCDATA / "DEEP_CHEM" / "toxcast.csv")
	sider = pd.read_csv(dirs.PROCDATA / "DEEP_CHEM" / "sider.csv")
	
	clintox = process_deepchem_data(colData,clintox,"Clinical Tox Result")
	sider = process_deepchem_data(colData,sider, "Side Effect")
	tox21 = process_deepchem_data(colData,tox21,"Tox Assay")
	toxcast = process_deepchem_data(colData,toxcast, "Tox Assay") 
	

	clintox.to_csv(dirs.PROCDATA / "experiments" / "clintox.csv", index = False)
	sider.to_csv(dirs.PROCDATA / "experiments" / "sider.csv",index = False)
	tox21.to_csv(dirs.PROCDATA / "experiments" / "tox21.csv",index=False)
	toxcast.to_csv(dirs.PROCDATA / "experiments" / "toxcast.csv",index=False)
if __name__ == '__main__':

	main()

	# parser = argparse.ArgumentParser(
	# 	prog = 'make_ColData ',
	# 	description = "Generate the colData and parse other annotationdb data")
	# parser.add_argument('-u', help = "route to all compounds in database")
	# parser.add_argument('-l', help = "lincs")
	# parser.add_argument('-j',help = 'jump cp')
	# parser.add_argument('-b',help = 'blood brain barrier')
	

	# args = parser.parse_args()

	# main(
	# 	db_url = args.u,
	# 	lincs_file = args.l,
	# 	jump_cp_file = args.j,
	# 	bbbp_file = args.b)
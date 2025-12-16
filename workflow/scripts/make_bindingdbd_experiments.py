from damply import dirs 
import pandas as pd
from collections import defaultdict
import numpy as np 
# def main(
# 	colData: pd.DataFrame,
# 	bdb_data: pd.DataFrame
# 	) -> None:
	


# def compute_affinity(ki,kd):
# 	if isinstance(ki,str):

def is_numeric(x):
	
	numeric_type = isinstance(x,int) or isinstance(x,float)
	return numeric_type

def main():
	colData = pd.read_csv(dirs.PROCDATA / "colData.csv",usecols = ['Pubchem CID', 'SMILES'])
	binding_db = pd.read_csv(dirs.RAWDATA / "BINDING_DB" / "BindingDB_All_202512_cleaned.csv")
	binding_db = binding_db.dropna(subset = ['Ki (nM)','Kd (nM)'], how= 'all')
	binding_db = binding_db[binding_db['PubChem CID'].isin(colData['Pubchem CID'])]
	targets = list(pd.unique(binding_db['Target Name']))
	cpds = list(pd.unique(binding_db['PubChem CID']))
	
	num_cpds = len(cpds)
	num_tgts = len(targets)
	
	
	cpd_to_idx = {cpds[i].item():i for i in range(num_cpds)}
	target_to_idx = {targets[i]:i for i in range(num_tgts)}
	seen_targets = []
	res = defaultdict(list)
	for cpd in cpds: 
		cpd_data = num_tgts*[np.nan]
		# print(len(cpd_data))


		temp = binding_db[binding_db['PubChem CID']==cpd]
		for idx, row in temp.iterrows():
			seen_targets.append(row['Target Name'])
			aff_col = 'Ki (nM)' if is_numeric(row['Kd (nM)']) else 'Kd (nM)'
			aff = row[aff_col]
			if isinstance(aff,str):
				if ">" in aff or "<" in aff:
					continue
				else:
					aff = float(aff)
			cpd_data[target_to_idx[row['Target Name']]]=aff
		res[cpd] = cpd_data
	# print(set(seen_targets))
	res = pd.DataFrame(res,index = targets)
	res = res.dropna(axis = 1,how='all')
	res = res.dropna(axis = 0,how='all')
	res.to_csv(dirs.PROCDATA / "experiments" / "binding_db.csv")
			
if __name__ == '__main__':
	main()
# if __name__ == '__main__':
# 	parser = argparse.ArgumentParser(
# 		prog = 'make_bdb_experiments ',
# 		description = "Generate the colData and parse other annotationdb data")
# 	parser.add_argument('-c', help = "coldata")
# 	parser.add_argument('-b', help = "bdbd")

# 	args = parser.parse_args()

# 	main(
# 		colData = args.c,
# 		bdb_data = args.b
# 		)
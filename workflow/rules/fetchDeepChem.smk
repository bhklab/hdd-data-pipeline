import pandas as pd
from pathlib import Path
from damply import dirs

deep_chem_urls = config['deep_chem']['urls']
deepchem_subdir = config['deep_chem']['subdir']
		


rule download_DeepChem:	
	output:
		bbbp = dirs.PROCDATA /deepchem_subdir/ "blood_brain_barrier.csv",
		toxcast = dirs.PROCDATA / deepchem_subdir / "toxcast.csv",
		sider = dirs.PROCDATA / deepchem_subdir / "sider.csv",
		clintox = dirs.PROCDATA / deepchem_subdir / "clintox.csv",
		tox21 = dirs.PROCDATA / deepchem_subdir  / "tox21.csv"

	run: 
		outpath = dirs.PROCDATA / deepchem_subdir
		Path(outpath).mkdir(parents=True,exist_ok=True)
		
		for dataset in deep_chem_urls.keys():
			#print(dataset)
			url = deep_chem_urls[dataset]
			
			if url[-3:]==".gz":
				compression = "gzip"
			else:
				compression = 'infer'

			data = pd.read_csv(url, compression=compression)
			data.to_csv(outpath / f"{dataset}.csv",index=False)
			





    	

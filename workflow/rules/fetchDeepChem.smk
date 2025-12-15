import pandas as pd
from damply import dirs

deep_chem_urls = config['deep_chem']['urls']
subdir = config['deep_chem']['subdir']
		


rule download_DeepChem:	
	output:
		bbbp = dirs.PROCDATA / subdir/ "blood_brain_barrier.csv",
		toxcast = dirs.PROCDATA / subdir / "toxcast.csv",
		sider = dirs.PROCDATA / subdir / "sider.csv",
		clintox = dirs.PROCDATA / subdir / "clintox.csv",
		tox21 = dirs.PROCDATA / subdir  / "tox21.csv"

	run: 
		
		outpath = dirs.PROCDATA / subdir
		Path(outpath).mkdir(parents=True,exist_ok=True)
		
		for dataset in  params.deep_chem_urls.keys():
			
			url = deep_chem_urls[dataset]
			
			if url[-3:]==".gz":
				compression = "gzip"
			else:
				compression = 'infer'

			data = pd.read_csv(url, compression=compression)
			data.to_csv(outpath / f"{dataset}.csv",index=False)
			





    	

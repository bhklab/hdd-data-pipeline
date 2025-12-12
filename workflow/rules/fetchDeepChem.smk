import pandas as pd
from damply import dirs





rule download_DeepChem:
	input:
		deep_chem_urls = config['deep_chem']['urls'],
		subdir_name = "DEEP_CHEM"

	output:
		bbbp = dirs.PROCDATA / input.subdir_name/ "blood_brain_barrier.csv"
		toxcast = dirs.PROCDATA / input.subdir_name/ "toxcast.csv"
		sider = dirs.PROCDATA / input.subdir_name/ "sider.csv"
		clintox = dirs.PROCDATA / input.subdir_name/ "clintox.csv"

	run: 
		outpath = dirs.PROCDATA / input.subdir_name
		Path(outpath).mkdir(parents=True,exist_ok=True)
		for dataset in  deep_chem_urls.keys():
			url = deep_chem_urls[dataset]
			if url[-3:]==".gz":
				compression = "gzip"
			else:
				compression = 'infer'

			data = pd.read_csv(url, compression=compression)
			data.to_csv(outpath / f"{dataset}.csv",index=False)
			





    	
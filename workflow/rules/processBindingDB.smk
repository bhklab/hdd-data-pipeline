from damply import dirs
import pandas as pd
import requests

info = config['binding_db']
subset = info['subset']
version = info['version']

rule download_BindingDB:
	input:
		base_url = info['base']

	output:
		zipped_data = dirs.RAWDATA / "BINDING_DB" / f"BindingDB_{subset}_{version}_tsv.zip"

	run: 
		outpath = dirs.RAWDATA / "BINDING_DB"
		fname = f"BindingDB_{subset}_{version}_tsv.zip"
		
		url = "_".join([input.base_url,subset, version,"tsv.zip"])

		Path(outpath).mkdir(parents=True,exist_ok=True)
	
		res = requests.get(url)
		
		with open(outpath / fname,"wb") as f
			f.write(res.content)


rule process_BindingDB:
	input: 
		raw_zip = dirs.RAWDATA / "BINDING_DB" / f"BindingDB_{subset}_{version}_tsv.zip",
		organisms = info['processing']['keep_organisms'],
		org_col = info['processing']['organism_col'],
		useful_cols = info['processing']['keep_cols'],
		assay_col = info['processing']['assay_col'],
		cid_col = info['processing']['cid_col'],
		output_cols = info['processing']['output_cols']
	

	output:
		cleaned_data = dirs.RAWDATA / "BINDING_DB" / f"BindingDB_{subset}_{version}_cleaned.csv"

		
	run: 		
		data = pd.read_csv(input.raw_zip,compression = 'zip',usecols = input.useful_cols)
	
		data = data[(data[input.processing.org_col].isin(input.processing.organisms)) & (data[input.assay_col].isnull())][input.output_cols]
	
		data.dropna(inplace=True, subset =[input.cid_col])
		
		data[input.cid_col] =[int(cid) for cid in data[input.cid_col].values]
		
		data.to_csv(output.cleaned_data,index=False)




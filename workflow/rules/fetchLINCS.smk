import pandas as pd
from pathlib import Path
from damply import dirs


lincs_version = config['lincs']['version']
lincs_subdir = config['lincs']['subdir'] 

rule download_LINCS:
	params:
		download_url = config['lincs']['url']
	
	output:
		lincs_raw = dirs.RAWDATA /lincs_subdir/ lincs_version / "compounds_raw.csv"
		
	run:
		Path(dirs.RAWDATA / lincs_subdir / lincs_version).mkdir(parents=True,exist_ok=True)
		compound_data = pd.read_csv(params.download_url, sep = "\t")
		compound_data.to_csv(output.lincs_raw,index=False)

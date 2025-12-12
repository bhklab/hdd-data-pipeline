import pandas as pd
from pathlib import Path
from damply import dirs


version = config['lincs']['version']


rule download_LINCS:
	input:
		download_url = config['lincs']['url']
	
	output:
		lincs_raw = dirs.RAWDATA / version / "compounds_raw.csv"
		
	run:
		Path(dirs.RAWDATA / version).mkdir(parents=True,exist_ok=True)
		compound_data = pd.read_csv(input.download_url, sep = "\t")
		compound_data.to_csv(output.lincs_raw,index=False)
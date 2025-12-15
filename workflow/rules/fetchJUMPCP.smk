from damply import dirs
import pandas as pd
import requests


subdir = config['jump_cp']['subdir']
version = config['jump_cp']['version']

rule download_JUMPCP:
	params:
		cpd_url = config['jump_cp']['url'],
	output: 
		data = dirs.RAWDATA / subdir  / version  /"JUMP_CP_compounds.csv"

	run:
		outpath = dirs.RAWDATA / subdir / version
		Path(outpath).mkdir(parents=True,exist_ok=True)
		

		data = pd.read_csv(params.cpd_url,compression='gzip')
		data.to_csv(outpath /"JUMP_CP_compounds.csv",index = False )

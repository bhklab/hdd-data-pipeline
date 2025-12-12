from damply import dirs
import pandas as pd
import requests




rule download_JUMPCP:
	input:
		cpd_url = config['jump_cp']['url'],
		subdir: "JUMP_CP",
		version = config['jump_cp']['version']

	output: 
		data = dirs.RAWDATA / input.subdir/ input.version/"JUMP_CP_compounds.csv"

	run:
		outpath = dirs.RAWDATA / input.subdir /input.version
		Path(outpath).mkdir(parents=True,exist_ok=True)
		

		data = pd.read_csv(cpd_url,compression='gzip')
		data.to_csv(outpath /"JUMP_CP_compounds.csv",index = False )
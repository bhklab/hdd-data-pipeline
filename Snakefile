import pandas as pd
import requests
from pathlib import Path
from damply import dirs


LINCS_VERSION = "LINCS_2020"
JUMP_CP_VERSION = "cpg0016"
BINDING_DB_VERSION = "202512"



rule all:
	input:
		dirs.RAWDATA


rule download_lincs:
	output:
		lincs_raw = dirs.RAWDATA / LINCS_VERSION / "compounds_raw.csv"
	
	run:
		import pandas as pd
		
		Path(dirs.RAWDATA / LINCS_VERSION).mkdir(parents=True,exist_ok=True)
		
		compound_url = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt"
		compound_data = pd.read_csv(url, sep = "\t")
		compound_data.to_csv(output.lincs_raw,index=False)



rule download_deepchem:
	output:
		output_file = dirs.RAWDATA / "DEEP_CHEM"  / "BBBP.csv" 

	run:
		import pandas as pd
		import gzip
		
		Path(dirs.RAWDATA / "DEEP_CHEM").mkdir(parents=True,exist_ok=True)
		bbbp_url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv"

		bbbp_data = pd.read_csv(bbbp_url)
		bbbp_data.to_csv(output_file,index=False)
		

		# saved for future use
		#toxcast_url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/toxcast_data.csv.gz"
		#sider_url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/sider.csv.gz",
		#clintox_url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/clintox.csv.gz"


rule download_jump_cp:
	output:
		jump_cp_file  = dirs.RAWDATA / "JUMP_CP" / JUMP_CP_VERSION / "compounds.csv"

	run: 
		import pandas as pd
		

		outpath = dirs.RAWDATA / "JUMP_CP" / JUMP_CP_VERSION
	
		Path(outpath).mkdir(parents=True,exist_ok=True)
		

		jcp_cpd_url = "https://github.com/jump-cellpainting/datasets/raw/refs/heads/main/metadata/compound.csv.gz"
		jcp_data = pd.read_csv(jcp_cpd_url,compression='gzip')
		jcp_data.to_csv(outpath /"compounds.csv", index = False)




rule download_binding_db:
	output:
		dirs.RAWDATA / "BINDING_DB"/ BINDING_DB_VERSION / "binding_db.csv"


	run:
		import pandas as pd
		import zipfile
		import os

		outpath = dirs.RAWDATA / "BINDING_DB"/ BINDING_DB_VERSION 
		Path(outpath).mkdir(parents=True,exist_ok=True)
		

		url = f"https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.jsp?download_file=/rwd/bind/downloads/BindingDB_All_{BINDING_DB_VERSION}_tsv.zip"
		res = requests.get(url)
		
		with open(outpath / "binding_db.csv.zip","wb") as f
			f.write(res.content)

		with zipfile.ZipFile(outpath /  "binding_db.csv.zip", 'r') as zf:
			zf.extract(outpath)

		os.system(f"rm {outpath / "binding_db.csv.zip"}")



rule process_data:
	input: 



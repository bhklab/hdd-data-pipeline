import sys
from pathlib import Path 
import pandas as pd
import urllib.request
import os
from typing import Union, List, Tuple
from collections import namedtuple
from damply import dirs
import json 
import tqdm 
import pubchempy as pcp

def fetch_LINCS(
	)->pd.DataFrame:
	"""
	
	"""
	compound_url = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/compoundinfo_beta.txt"
	# sig_url = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt"

	write_dir = dirs.RAWDATA / "LINCS"
	write_dir.mkdir(exist_ok=True)
	
	for url in [compound_url]:
		output_file = url.split("/")[-1]
		output_path  = write_dir /output_file
		
		if output_file not in list([x.name for x in write_dir.glob("*")]):
			os.system(f"wget {url} -P {str(write_dir)}")



	
	cpd_data = pd.read_csv(write_dir / "compoundinfo_beta.txt",sep="\t", usecols = ['pert_id', 'cmap_name','canonical_smiles','inchi_key'])
	cpd_data.dropna(inplace=True)
	cpd_data.drop_duplicates(inplace=True)
	# cpd_data = cpd_data.iloc[:10,]
	cids, smiles = [],[]
	keep_idx = []
	for idx, row in tqdm.tqdm(cpd_data.iterrows(),total = cpd_data.shape[0],leave=False):
		try:
			cpd = pcp.get_compounds(row['canonical_smiles'],'smiles')[0]
			cids.append(cpd.cid)
			smiles.append(cpd.canonical_smiles)
			keep_idx.append(idx)
		except:
			continue
	cpd_data = cpd_data.iloc[keep_idx,:]
	cpd_data['fetched_smiles'] = smiles
	cpd_data['PubChem CID'] = cids
	cpd_data.to_csv(dirs.PROCDATA / "lincs_annotated_raw.csv")


	
def fetch_DeepChem()->pd.DataFrame:
	urls = ["https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv",
		"https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz",
		"https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/toxcast_data.csv.gz",
		"https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/sider.csv.gz",
		"https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/clintox.csv.gz"]
	# pass#"https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv"

def fetch_JUMPCP()->pd.DataFrame:
	url = ""
	



def fetch_binding_db():
	url = "https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202511_tsv.zip"
	

def fetch_drug_universe():
	fetch_LINCS()




if __name__ == '__main__':
	fetch_drug_universe()
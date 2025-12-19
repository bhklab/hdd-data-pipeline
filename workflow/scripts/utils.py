from collections import defaultdict
from typing import Dict, List, Union

import pandas as pd

GOLD_STANDARD_AIDS = [
	485290,
	1508612,
	1645840,
	1645841,
	1645842,
	492947,
	1030,
	743075,
	743080,
	588795,
	2101,
	602179,
	504327,
	995,
	493208,
	1777,
	1631,
	743094,
	651631,
	504847,
	1159551,
	1259242,
	1259241,
	743012,
	1224868,
	1224880,
	1346977,
	743014,
	1224870,
	1224872,
	1224874,
	1224877,
	1224885,
	1224886,
	743015,
	1224873,
	1224887,
	1224889,
	1224867,
	720516,
	651632,
	651634,
]


def process_single_drug(
	drug_info,
	drug_details: Dict[str, Union[int, float, str]],
	colData: defaultdict(list),
	all_bioassays: Dict[str, Dict],
	seen_bioassays: List[int],
	lincs_compounds: pd.DataFrame,
	jump_cp_compounds: pd.DataFrame,
	blood_brain_perm: pd.DataFrame,
	cids: List,
) -> None:
	# Molecule Name

	colData['Molecule Name'].append(drug_info['name'])
	colData['Pubchem CID'].append(drug_info['cid'])
	colData['InChIKey'].append(drug_info['inchikey'])
	colData['SMILES'].append(drug_info['smiles'])
	cids.append(drug_info['cid'])
	# store for later use
	smiles_str = drug_info['smiles']
	cid = drug_info['cid']

	colData['Molecular Formula'].append(drug_details['molecular_formula'])
	colData['IUPAC Name'].append(drug_details['iupac_name'])
	colData['ChEMBL ID'].append(drug_details['molecule_chembl_id'])
	colData['PubChem 2D Fingerprint'].append(drug_details['fingerprint_2d'])

	# MOA and Approval
	mechanisms = drug_details['mechanisms']

	if len(mechanisms) == 0:
		colData['Mechanism of Action'].append('None')
	else:
		colData['Mechanism of Action'].append(
			drug_details['mechanisms'][0]['mechanism_of_action']
		)

	colData['FDA Approved'].append(drug_details['fda_approval'])

	## Add in Molecular Information (molecular weight + Lipinski Filters)
	## 	for the curious: https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five

	colData['Molecular Weight'].append(drug_details['molecular_weight'])
	colData['XlogP'].append(drug_details['xlogp'])
	colData['Hydrogen Bond Donors'].append(drug_details['h_bond_donor_count'])
	colData['Hydrogen Bond Acceptors'].append(drug_details['h_bond_acceptor_count'])
	colData['Exact Molecular Mass'].append(drug_details['exact_mass'])
	colData['DILI Severity'].append(drug_details['toxicity']['dili_severity_grade'])
	colData['DILI Annotation'].append(drug_details['toxicity']['dili_annotation'])
	colData['Hepatotoxicity Likelihood (Detailed)'].append(
		drug_details['toxicity']['hepatotoxicity_likelihood_score']
	)
	score = (
		drug_details['toxicity']['hepatotoxicity_likelihood_score']
		.split(':')[1]
		.lstrip()
		.split()[0]
	)
	colData['Hepatotoxiciy Likelihood (Score)'].append(score)

	## Check Against The Broad Data

	# print("pingo herebo")
	# print(lincs_compounds)
	# print(lincs_compounds['inchi_key'].value_counts())
	l1k_subset = lincs_compounds[lincs_compounds['inchi_key'] == drug_info['inchikey']]
	jump_subset = jump_cp_compounds[
		jump_cp_compounds['Metadata_InChIKey'] == drug_info['inchikey']
	]

	if l1k_subset.shape[0] == 0:
		# print("thrig  plibbus")
		colData['In L1000'].append(False)
		colData['L1000 ID'].append('-')
	else:
		colData['In L1000'].append(True)
		colData['L1000 ID'].append(l1k_subset['pert_id'].values[0])

	if jump_subset.shape[0] == 0:
		colData['In JUMP-CP'].append(False)
		colData['JUMP-CP ID'].append('-')
	else:
		colData['In JUMP-CP'].append(True)
		colData['JUMP-CP ID'].append(jump_subset['Metadata_JCP2022'].values[0])

	blood_brain = blood_brain_perm[blood_brain_perm['smiles'] == smiles_str]
	if blood_brain.shape[0] == 0:
		colData['BBB Permeable'].append('Unknown')
	else:
		colData['BBB Permeable'].append(blood_brain['p_np'].values[0])

	# Cache bioassays for post-processing
	all_bioassays[drug_details['cid']] = drug_details['bioassays']
	seen_bioassays.extend([assay['aid'] for assay in drug_details['bioassays']])






def main():


if __name__ == '__main__':
	parser.add_argument('-r', help = "radius list")
	parser.add_argument('-d', help = "dimension list")main()


	radius_list = args.r,
		dim_list = args.d,	
	# make the fingerprint generators
	fingerprint_generators = make_fingerprint_generators(radius_list, dim_list)
	
	# store as a dictionary for breaking out as files later
	fp_data = {k:defaultdict(list) for k in fingerprint_generators}

	fp_data = fp_data,
	fingerprint_generators = fingerprint_generators
	fp_data:Dict[str, defaultdict(list)],

	fingerprint_generators: Dict[str,rdkit.Chem.rdFingerprintGenerator]

	fingerprint_generators: Dict[str,rdkit.Chem.rdFingerprintGenerator]

	# Make the fingerpints
	rdk_mol = Chem.MolFromSmiles(smiles_str)
	for fp_type in fingerprint_generators:
		generator = fingerprint_generators[fp_type]
		fp = generator.GetCountFingerprintAsNumPy(rdk_mol)
		fp_data[fp_type]['CID'].append(cid)
		for dim in range(fp.shape[0]):
			fp_data[fp_type][f'V{dim+1}'].append(fp[dim].item())
		
	
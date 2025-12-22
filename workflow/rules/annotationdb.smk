from damply import dirs


rule fetch_from_AnnotationDB:
	params:
		db_url = config["colData"]["db_url"]

	input:
		raw_data = rules.fetch_AnnotationDB_raw.output.raw,
		lincs_file = rules.download_LINCS.output.lincs_raw,
		jump_file = rules.download_JUMPCP.output.data,
		bbbp_file = rules.download_DeepChem.output.bbbp

	output:
		colData = dirs.PROCDATA / "colData.csv",
		bioassays = dirs.PROCDATA / "experiments" / "bioassays.csv"

	threads: 1

	shell:
		"""
		mkdir -p {dirs.PROCDATA} {dirs.PROCDATA}/experiments
		PYTHONPATH=workflow/scripts python ./workflow/scripts/process_annotationdb.py \
			-i {input.raw_data} -l {input.lincs_file} -j {input.jump_file} -b {input.bbbp_file}
		"""


rule fetch_AnnotationDB_raw:
	params:
		db_url = config["colData"]["db_url"]

	output:
		raw = dirs.RAWDATA / "ANNOTATION_DB" / "compound_details.jsonl"

	threads: 1

	shell:
		"""
		mkdir -p {dirs.RAWDATA}/ANNOTATION_DB
		PYTHONPATH=workflow/scripts python ./workflow/scripts/fetch_annotationdb.py -u {params.db_url} -o {output.raw}
		"""

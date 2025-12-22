from damply import dirs


rule construct_MAE:
	input:
		colData = rules.fetch_from_AnnotationDB.output.colData,
		bioassays = rules.fetch_from_AnnotationDB.output.bioassays,
		binding_db = rules.make_bindingdb_experiments.output.binding_db,
		deepchem = rules.make_deepchem_experiments.output.experiments,
		fingerprints = rules.make_fingerprints.output.fingerprints

	output:
		mae = dirs.RESULTS / "hdd.RDS"

	shell:
		"""
		mkdir -p {dirs.RESULTS}
		Rscript ./workflow/scripts/construct_MAE.R
		"""

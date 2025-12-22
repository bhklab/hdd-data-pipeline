from damply import dirs


rule make_bindingdb_experiments:
	input:
		colData = rules.fetch_from_AnnotationDB.output.colData,
		bdb_data = rules.process_BindingDB.output.cleaned_data

	output:
		binding_db = dirs.PROCDATA / "experiments" / "binding_db.csv"

	shell:
		"""
		mkdir -p {dirs.PROCDATA}/experiments
		python ./workflow/scripts/make_bindingdbd_experiments.py -c {input.colData} -b {input.bdb_data}
		"""

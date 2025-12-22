from damply import dirs

deepchem_subdir = config["deep_chem"]["subdir"]
deepchem_experiments = ["toxcast", "tox21", "sider", "clintox"]


rule make_deepchem_experiments:
	params:
		deepchem_subdir = deepchem_subdir

	input:
		colData = rules.fetch_from_AnnotationDB.output.colData,
		deepchem_files = expand(
			dirs.PROCDATA / deepchem_subdir / "{dataset}.csv",
			dataset=deepchem_experiments,
		)

	output:
		experiments = expand(
			dirs.PROCDATA / "experiments" / "{dataset}.csv",
			dataset=deepchem_experiments,
		)

	shell:
		"""
		mkdir -p {dirs.PROCDATA}/experiments
		python ./workflow/scripts/make_deepchem_experiments.py -c {input.colData} -s {params.deepchem_subdir}
		"""

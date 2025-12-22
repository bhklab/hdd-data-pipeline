from pathlib import Path
from damply import dirs


configfile: "config/pipeline.yaml"

### these are the downloading processes
include: "workflow/rules/processBindingDB.smk"
include: "workflow/rules/fetchLINCS.smk"
include: "workflow/rules/fetchJUMPCP.smk"
include: "workflow/rules/fetchDeepChem.smk"

### rules for processing the assay files


## this will also pull and make the bioassay experiments

rule fetch_from_AnnotationDB:
	params:
		db_url = config['colData']['db_url']

	input:
		raw_data = rules.fetch_AnnotationDB_raw.output.raw,
		lincs_file = rules.download_LINCS.output.lincs_raw,
		jump_file = rules.download_JUMPCP.output.data,
		bbbp_file = rules.download_DeepChem.output.bbbp,

	output:
		colData = dirs.PROCDATA/ "colData.csv",
		bioassays = dirs.PROCDATA / "experiments" /  "bioassays.csv"

	threads: 1

	shell:
		"""
		mkdir -p {dirs.PROCDATA} {dirs.PROCDATA}/experiments
		PYTHONPATH=workflow/scripts python ./workflow/scripts/process_annotationdb.py \
			-i {input.raw_data} -l {input.lincs_file} -j {input.jump_file} -b {input.bbbp_file}
		"""


rule fetch_AnnotationDB_raw:
	params:
		db_url = config['colData']['db_url']

	output:
		raw = dirs.RAWDATA / "ANNOTATION_DB" / "compound_details.jsonl"

	threads: 1

	shell:
		"""
		mkdir -p {dirs.RAWDATA}/ANNOTATION_DB
		PYTHONPATH=workflow/scripts python ./workflow/scripts/fetch_annotationdb.py -u {params.db_url} -o {output.raw}
		"""
		
rule all:
	input: 
		rules.fetch_from_AnnotationDB.output.colData,
		rules.fetch_from_AnnotationDB.output.bioassays
		
# rule make_fingerprints:
# 	params:
# 		radius_list = config['colData']['fingerprints']['radius_list'],
# 		dim_list = config['colData']['fingerprints']['dim_list']

# 	input:
# 		rules.fetch_from_AnnotationDB.output.colData

# 	output:
# 		expand(dirs.PROCDATA / "experiments"/ "Morgan({rad},{dim}).csv",rad = params.radius_list, dim = params.dim_list)
	
# 	shell:
# 		"python3 ./workflow/scripts/make_fingerprints.py -r {params.radius_list} -d {params.dim_list}"



# rule make_bindingdb_experiment:
# 	input: 
# 		colData = rules.fetch_from_AnnotationDB.output.colData,
# 		bdb_data = rules.process_BindingDB.output.cleaned_data

# 	output:
# 		bdb_experiment = dirs.PROCDATA / "experiments"/ "binding_db.csv"


# 	shell:
# 		"python3 ./workflow/scripts/make_bindingdb_experiment.py -c {input.colData} -b {input.bdb_data}"


# rule make_deepchem_experiments:
# 	params: 
# 		exps = ['toxcast', 'tox21', 'sider','clintox'],
# 		subdir = config['deep_chem']['subdir']
	
# 	input:
# 		colData = rules.fetch_from_AnnotationDB.output.colData,
# 		expand(dirs.PROCDATA/ subdir/"{dataset}.csv",dataset = params.exps)
	
# 	output: 
# 		expand(dirs.PROCDATA / "experiments" /"{dataset}.csv",dataset = params.exps)



# rule make_multi_assay_experiment:
# 	input:
# 		rules.fetch_from_AnnotationDB.output.colData,
# 		rules.fetch_from_AnnotationDB.output.bioassays,
# 		rules.make_fingerprints.output

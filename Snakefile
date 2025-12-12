from pathlib import Path
from damply import dirs


configfile: "config/pipeline.yaml"

### these are the downloading processes
include: "workflow/rules/processBindingDB.smk"
include: "workflow/rules/fetchLINCS.smk"
include: "workflow/rules/fetchJUMPCP.smk"
include: "workflow/rules/fetchDeepChem.smk"

### rules for processing the assay files

rule all:
	input: 
		dirs.PROCDATA / "colData.csv",
		dirs.PROCDATA /"bioassays.csv",
		rules.download_LINCS.output.lincs_raw,
		rules.download_JUMPCP.output.data,
		rules.process_BindingDB.output.cleaned_data



rule fetch_from_AnnotationDB:
	params:
		radius_list = config['colData']['fingerprints']['radius_list'],
		dim_list = config['colData']['fingerprints']['dim_list'],
		db_url = config['colData']['db_url']
	

	input:
		lincs_file = rules.download_LINCS.output.lincs_raw,
		jump_file = rules.download_JUMPCP.output.data,

	
	output:
		colData = dirs.PROCDATA/ "colData.csv",
		bioassays = dirs.PROCDATA / "bioassays.csv"

	shell:
		"python3 ./workflow/scripts/make_colData.py -u {params.db_url} -r {params.radius_list} -d {params.dim_list} -l {input.lincs_file} -j {input.jump_file}"
		




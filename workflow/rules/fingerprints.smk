from damply import dirs

fingerprint_radii = config["colData"]["fingerprints"]["radius_list"]
fingerprint_dims = config["colData"]["fingerprints"]["dim_list"]


rule make_fingerprints:
	params:
		radius_list = repr(fingerprint_radii),
		dim_list = repr(fingerprint_dims)

	input:
		rules.fetch_from_AnnotationDB.output.colData

	output:
		fingerprints = expand(
			dirs.PROCDATA / "experiments" / "fingerprints" / "Morgan.{rad}.{dim}.csv",
			rad=fingerprint_radii,
			dim=fingerprint_dims,
		)

	shell:
		"""
		PYTHONPATH=workflow/scripts python -c "import make_fingerprints as mf; mf.main(radius_list={params.radius_list}, dimension_list={params.dim_list})"
		"""

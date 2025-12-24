# Changelog

## [0.2.0](https://github.com/bhklab/hdd-data-pipeline/compare/v0.1.0...v0.2.0) (2025-12-24)


### Features

* add basic plots to qc ([c6e4b30](https://github.com/bhklab/hdd-data-pipeline/commit/c6e4b3026bd974c5a514791c6befa799c01ba312))
* add qc ([6678a20](https://github.com/bhklab/hdd-data-pipeline/commit/6678a20ce4684b3e45403e956cdc9381f482bc13))
* batch annotationdb fetch ([3a67a15](https://github.com/bhklab/hdd-data-pipeline/commit/3a67a15c53137542c954a4d959e0f33622ee2e6e))
* export MAE to csv ([30fb249](https://github.com/bhklab/hdd-data-pipeline/commit/30fb2499926e2b475715f8415bb1123be75360b9))
* parameterize deepchem experiments inputs ([3911963](https://github.com/bhklab/hdd-data-pipeline/commit/3911963dec9922193fe7cfe917886413aab5278d))


### Bug Fixes

* add retries for annotationdb requests ([3cc3e20](https://github.com/bhklab/hdd-data-pipeline/commit/3cc3e2006c0f8015c837778bb0c74788701e4d41))
* address dropped compounds ([3621740](https://github.com/bhklab/hdd-data-pipeline/commit/362174095a6a7eb88028a99401b8581b8656fc49))
* allow NA in tox21 integer conversion ([fadcd51](https://github.com/bhklab/hdd-data-pipeline/commit/fadcd512ef5ed9eeca4f2e098f9d5e0a5cbce4c3))
* coerce coldata na values ([585f189](https://github.com/bhklab/hdd-data-pipeline/commit/585f189f1bc6f45db58e8683c02826582a5f4030))
* define annotationdb raw rule before use ([5fe4eed](https://github.com/bhklab/hdd-data-pipeline/commit/5fe4eed3c7cbd1527f5639b8a83c1e457ce58137))
* import Path in download rules ([91bf27b](https://github.com/bhklab/hdd-data-pipeline/commit/91bf27b33512c419ee760dad457682f564fb401d))
* make_colData atomic per compound to avoid mismatched col lengths ([6ea36f6](https://github.com/bhklab/hdd-data-pipeline/commit/6ea36f63bc398eb0903ee6271aa7f02bcf8e6df2))
* NA bioassays ([798c287](https://github.com/bhklab/hdd-data-pipeline/commit/798c287bb1c7844bb1b5298b2769701136f97216))
* parameterize bindingdb experiment inputs ([3c875dc](https://github.com/bhklab/hdd-data-pipeline/commit/3c875dc52208eada660e683f3f798d0b985cb5ee))
* reduce annotationdb batch size ([8f1d756](https://github.com/bhklab/hdd-data-pipeline/commit/8f1d756fcfcbd2b3b929214a5a00c840621d238b))
* skip filtered assays ([c69c867](https://github.com/bhklab/hdd-data-pipeline/commit/c69c867c531916a280787ca66b687e57f3ab6a47))
* split annotationdb batches on server errors ([a424d09](https://github.com/bhklab/hdd-data-pipeline/commit/a424d09a431d76e4560bef2f473609f56d74c01a))

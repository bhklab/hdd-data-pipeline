library(MultiAssayExperiment)


exp.names <- c("BindingDB",'Bioassays', 'ToxCast', 'Tox21', 'ClinTox','Sider')
colData <- read.csv("data/procdata/colData.csv")
# fps <- read.csv("data/procdata/experiments/fingerprints.csv")
bdb <- read.csv("data/procdata/experiments/binding_db.csv",check.names=FALSE)
bioassays <- read.csv("data/procdata/experiments/filtered_bioassays.csv",check.names=FALSE)
toxcast <- read.csv("data/procdata/experiments/toxcast.csv",check.names = FALSE)
tox21 <- read.csv("data/procdata/experiments/tox21.csv",check.names = FALSE)
clintox <- read.csv("data/procdata/experiments/clintox.csv",check.names = FALSE)
sider <- read.csv("data/procdata/experiments/sider.csv",check.names = FALSE)

el = ExperimentList(
	 BindingDB= bdb,
	 Bioassays = bioassays,
	 ToxCast = toxcast,
	 Tox21 = tox21,
	 ClinTox = clintox,
	 SIDER = sider
	 )

sml <- lapply(el, function(se) { 
	data.frame(
    primary = colnames(se),
    colname = colnames(se),
    stringsAsFactors = FALSE
  )
})

mae = MultiAssayExperiment(experiments = el, colData = colData, sampleMap = listToMap(sml))
saveRDS(mae,"data/results/hdd.RDS")
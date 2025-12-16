import pandas as pd

from damply import dirs 

import utils 

df = pd.read_csv(dirs.PROCDATA / "experiments" / "bioassays.csv")
# print(df)
keep_assays = [f"AID_{a}" for a in utils.GOLD_STANDARD_AIDS]
# print(keep_assays)
df = df[df['Assay'].isin(keep_assays)]
df.to_csv(dirs.PROCDATA / "experiments" / "filtered_bioassays.csv")
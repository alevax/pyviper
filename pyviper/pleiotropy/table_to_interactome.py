from typing import Dict
import pandas as pd

# optional progress bar
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x, **kwargs):
        return x
# Function to Convert Interactome to Regulon (TabletoRegulon)
def TableToInteractome(net_table: pd.DataFrame) -> Dict[str, Dict[str, pd.Series]]:
    # First block of TableToInteractome: initialize my_reg with a blank
    # structure for each regulator that appears in net_table.
    # 1) blank_reg in R -> we'll create a fresh copy per TF in Python
    #    (do NOT reuse the same object across TFs)
    def _fresh_blank():
        # make sure each TF gets a space for tfmode and likelihood
        return {"tfmode": [], "likelihood": [], "_targets": []}

    # Purpose first loop: just to initialize each regulator with a blank regulon
    # At the end of this, every TF has an empty tfmode and likelihood.
    my_reg: Dict[str, Dict[str, list]] = {}

    # loop over each row with progress bar
    for _, row in tqdm(net_table.iterrows(), total=len(net_table), desc="Initializing regulons"):
        regulator = row["regulator"]
        # Overwrite every time (like in R)
        my_reg[regulator] = _fresh_blank()

    # Purpose second loop: filling in the values from the net_table for each target and regulator
    for _, row in tqdm(net_table.iterrows(), total=len(net_table), desc="Filling regulons"):
        r = row["regulator"]
        t = row["target"]
        w = row["likelihood"]
        m = row["mor"]

        # In Python we collect targets, then set them as index (names) at the end
        if r not in my_reg:
            my_reg[r] = _fresh_blank()
        my_reg[r]["_targets"].append(str(t))
        my_reg[r]["likelihood"].append(float(w))
        my_reg[r]["tfmode"].append(float(m))

    # --- Finalize: convert lists to named Series (index = target), like R named vectors ---
    for r, reg in my_reg.items():
        targets = reg["_targets"]
        reg["likelihood"] = pd.Series(reg["likelihood"], index=targets, dtype=float)
        reg["tfmode"]     = pd.Series(reg["tfmode"],     index=targets, dtype=float)
        del reg["_targets"]  # cleanup helper list

    return my_reg

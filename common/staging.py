"""AJCC tumour stage parsing shared by the survival modules."""
from __future__ import annotations

import pandas as pd

# "Stage IIIA" -> "III", "Stage IVB" -> "IV": substages fold into their stage group
STAGE_PATTERN = r"^STAGE\s+(IV|I{1,3})[A-C]?$"


def ajcc_stage_group(values: pd.Series) -> pd.Series:
    """Map AJCC stage labels to "I", "II", "III" or "IV".

    Missing values stay missing. A label that is present but does not parse is
    reported rather than dropped silently.
    """
    text = values.astype("string").str.strip().str.upper()
    group = text.str.extract(STAGE_PATTERN)[0]
    unparsed = text[text.notna() & group.isna()]
    if len(unparsed):
        print(f"  [stage] {len(unparsed)} labels not recognised: "
              f"{unparsed.value_counts().to_dict()}")
    return group

# Exported from the PSCM workbench (record 'compareBounds2Models', status: Drafting, updated 2026-10-05T16:00:00Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""Compare reaction bounds between two models.

MATLAB source: src/analysis/wholeBody/PSCMToolbox/compareBounds2Models.m, commit: TODO
Translation status: Draft, not parity-tested.

cobrapy used
------------
model.reactions (DictList): .has_id and .get_by_id; Reaction.lower_bound and
Reaction.upper_bound (cobra.core.reaction).

Output format
-------------
Same four outputs as MATLAB. Cell arrays become lists:
- r1_missing, r2_missing: lists of reaction ids
- r12_bounds, r12_bounds_mismatch: lists of
  [reaction id, description, value in model1, value in model2]
Rows are in the same order as MATLAB: reactions sorted by id (unique() sorts),
lower bound row before upper bound row.

Deviations from MATLAB
----------------------
1. MATLAB builds R1_missing / R2_missing by vertically concatenating char arrays
   ([R1_missing; R12{i}]). That errors as soon as two missing ids have different
   lengths. This translation returns lists, so it does not fail there.
2. If no bound is identical, MATLAB never assigns R12_bounds and errors on
   return. This translation returns an empty list.
3. MATLAB displays mismatch_lb and mismatch_ub without semicolons; this
   translation prints them in the same 'name = value' form.

Comparisons use ==, as MATLAB does: NaN never equals NaN, and inf equals inf.
"""
from typing import List, Tuple

from cobra import Model

Row = List[object]


def compare_bounds_2_models(model1: Model, model2: Model
                            ) -> Tuple[List[str], List[str], List[Row], List[Row]]:
    """Find reactions missing from either model and compare bounds of shared ones.

    MATLAB: [R1_missing, R2_missing, R12_bounds, R12_bounds_mismatch] =
                compareBounds2Models(model1, model2)

    Returns
    -------
    r1_missing : reactions absent from model1 (present in model2)
    r2_missing : reactions absent from model2 (present in model1)
    r12_bounds : shared reactions with an identical lower/upper bound, with values
    r12_bounds_mismatch : shared reactions whose lower/upper bound differs, with values
    """
    # MATLAB: R12 = unique([R1; R2])  (sorted, duplicates removed)
    r12 = sorted({r.id for r in model1.reactions} | {r.id for r in model2.reactions})

    r1_missing: List[str] = []
    r2_missing: List[str] = []
    r12_bounds: List[Row] = []
    r12_bounds_mismatch: List[Row] = []
    mismatch_lb = 0
    mismatch_ub = 0

    for rid in r12:
        if not model1.reactions.has_id(rid):        # reaction does not exist in model1
            r1_missing.append(rid)
        elif not model2.reactions.has_id(rid):      # reaction does not exist in model2
            r2_missing.append(rid)
        else:
            rxn1 = model1.reactions.get_by_id(rid)
            rxn2 = model2.reactions.get_by_id(rid)
            # lower bound
            if rxn1.lower_bound == rxn2.lower_bound:
                r12_bounds.append([rid, "identical lower bound",
                                   rxn1.lower_bound, rxn2.lower_bound])
            else:
                r12_bounds_mismatch.append([rid, "NOT identical lower bound",
                                            rxn1.lower_bound, rxn2.lower_bound])
                mismatch_lb += 1
            # upper bound
            if rxn1.upper_bound == rxn2.upper_bound:
                r12_bounds.append([rid, "identical upper bound",
                                   rxn1.upper_bound, rxn2.upper_bound])
            else:
                r12_bounds_mismatch.append([rid, "NOT identical upper bound",
                                            rxn1.upper_bound, rxn2.upper_bound])
                mismatch_ub += 1

    print(f"mismatch_lb = {mismatch_lb}")
    print(f"mismatch_ub = {mismatch_ub}")
    return r1_missing, r2_missing, r12_bounds, r12_bounds_mismatch

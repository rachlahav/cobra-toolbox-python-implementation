# Exported from the PSCM workbench (record 'adjustWholeBodyRxnCoeff', status: Drafting, updated 2026-10-05T16:00:00Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""Adjust organ coefficients in the whole-body maintenance reaction.

MATLAB source: src/analysis/wholeBody/PSCMToolbox/adjustWholeBodyRxnCoeff.m
               (Ines Thiele, 2012-2020), commit: TODO
Translation status: Draft, not parity-tested.

cobrapy used
------------
Reaction.add_metabolites(..., combine=False) (cobra.core.reaction): replaces the
coefficient of a metabolite already in the reaction, adds it if absent, and
removes it if the new coefficient is 0. Reverted on exit inside `with model:`.

Deviations from MATLAB
----------------------
1. If 'Whole_body_objective_rxn' is missing, MATLAB silently changes nothing
   (empty column index). This translation also changes nothing, but warns.
2. MATLAB errors when listCoeff is shorter than listOrgan (index out of range).
   This translation raises ValueError in the same case instead of letting zip()
   silently drop organs.
3. The model.A / model.S handling is not translated. cobrapy stores one
   stoichiometry per reaction, so there is nothing to copy back.

Source quirks
-------------
- MATLAB ends with model.S = model.A. For a model that already has an old-style
  model.A = [S; C], this makes S include the coupling rows. Not reproduced: in
  cobrapy only the metabolite coefficients change.
- Organ metabolites are found by prefix match on '<Organ>_biomass' (strmatch)
  intersected with a substring match on '_dummy_objective' (strfind), so
  '<Organ>_biomassX..._dummy_objective' would also match. Kept for parity.

Note on zero coefficients: setting a coefficient to 0 removes the metabolite from
the reaction in cobrapy. In MATLAB, assigning 0 into the sparse S also removes
the stored entry, so both give the same stoichiometry.
"""
import warnings
from typing import Sequence

from cobra import Model

WHOLE_BODY_RXN = "Whole_body_objective_rxn"


def adjust_whole_body_rxn_coeff(model: Model, list_organ: Sequence[str],
                                list_coeff: Sequence[float]) -> Model:
    """Set each organ's coefficient in Whole_body_objective_rxn to -fraction * 100.

    MATLAB: model = adjustWholeBodyRxnCoeff(model, listOrgan, listCoeff)

    Parameters
    ----------
    model : cobra.Model
        Whole-body model. Changed in place and returned.
    list_organ : sequence of str
        Organ names, e.g. ["Liver", "Muscle"].
    list_coeff : sequence of float
        Organ weight fractions (0-1), in the same order as list_organ. They are
        written as percentages, as in MATLAB.
    """
    # MATLAB: wholeBodyRxnID = find(ismember(model.rxns, wholeBodyRxn))
    if not model.reactions.has_id(WHOLE_BODY_RXN):
        warnings.warn(f"{WHOLE_BODY_RXN} not in model; nothing changed")
        return model
    whole_body_rxn = model.reactions.get_by_id(WHOLE_BODY_RXN)

    if len(list_coeff) < len(list_organ):
        raise ValueError("list_coeff is shorter than list_organ "
                         f"({len(list_coeff)} < {len(list_organ)})")

    # MATLAB: organID2 = find(~cellfun(@isempty, strfind(model.mets, '_dummy_objective')))
    # The same for every organ, so it is computed once.
    dummy_mets = [met for met in model.metabolites if "_dummy_objective" in met.id]

    for organ, coeff in zip(list_organ, list_coeff):
        # MATLAB: organID = intersect(strmatch(strcat(organ, '_biomass'), model.mets), organID2)
        organ_mets = [met for met in dummy_mets if met.id.startswith(f"{organ}_biomass")]
        if organ_mets:
            # MATLAB: model.A(organID, wholeBodyRxnID) = -listCoeff(i)*100
            # (fractions are stored as percentages in the whole-body objective)
            whole_body_rxn.add_metabolites({met: -coeff * 100 for met in organ_mets},
                                           combine=False)
    return model

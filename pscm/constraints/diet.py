# Exported from the PSCM workbench (record 'setDietConstraints', status: Drafting, updated 2026-10-05T05:07:57.646Z).
# Edit in the workbench, then re-run tools/export_workbench.py.
"""Set diet constraints on a whole-body model.

MATLAB source: src/analysis/wholeBody/PSCMToolbox/setConstraints/setDietConstraints.m
               (Ines Thiele, 2016-2019), commit: TODO
Translation status: Draft, not parity-tested.

What it does
------------
Closes every Diet_EX_ uptake, opens a small uptake (-0.1) for metabolites AGORA
microbes need that the diet lacks, opens uptake for compounds the diet database
is missing, then constrains each diet reaction to between 1.2x and 0.8x of the
diet value (mmol/day/person), with special rules for micronutrients, ions and
sulfate.

Conventions applied
-------------------
- Reaction ids are used exactly as in the MATLAB model (Diet_EX_glc_D[d]).
- strmatch(prefix, ...) -> str.startswith; ismember -> exact set membership.
- Final (lb, ub) pairs are computed first and assigned once with
  reaction.bounds, because cobrapy raises as soon as lb > ub.
- The model is changed in place and returned (lab decision pending).
- SetupInfo lives in model.notes["SetupInfo"].

Deviations from MATLAB
----------------------
1. An unknown diet name raises ValueError. MATLAB falls through and fails later
   at Diet(:,1).
2. A diet row whose final lb > ub (only possible with a negative diet value)
   raises ValueError. MATLAB stores it and the LP becomes infeasible.
3. Diet rows naming reactions not in the model are skipped as in MATLAB, but
   reported with one warning instead of silently.
4. Diet values may be numbers or strings. MATLAB str2num needs strings.
5. The unused copy `modelO = model` is not translated.

Source quirks kept for parity (see workbench notes)
---------------------------------------------------
- Micronutrients: values <= 10 get lb = -10*factor; the "> 0.1" branch is only
  reached for values > 10 and gives lb = -120*value*factor.
- The three upper-bound branches are identical.
- MissingDietCompounds and Diet_EX_chol[d] are set even when the diet contains
  them; the diet rows then overwrite the lower bound.
- Upper bounds are negative, so diet uptake is enforced (the MATLAB `if 1` block
  is commented "set if statement to 0" to not enforce it).
"""
import csv
import warnings
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple, Union

from cobra import Model

from pscm.data.agora_essential import AGORA_ESSENTIAL

DietTable = List[Tuple[str, str]]

DIET_DIR = Path(__file__).resolve().parent.parent / "data" / "diets"

# Diet names accepted by the MATLAB function, mapped to the script each one runs.
DIET_FILES = {
    "EUAverageDiet": "EUAverageDietNew",
    "HighFiberDiet": "HighFiberDiet",
    "HighProteinDiet": "HighProteinDiet",
    "UnhealthyDiet": "UnhealthyDiet",
    "VegetarianDiet": "VegetarianDiet",
}

# MATLAB: microbiotaEnabling = 1 (hard-coded)
MICROBIOTA_ENABLING = True

# MATLAB: MissingDietCompounds, first list (lb = -50)
MISSING_DIET_COMPOUNDS_1 = [
    "Diet_EX_asn_L[d]", "Diet_EX_gln_L[d]", "Diet_EX_chol[d]", "Diet_EX_crn[d]",
    "Diet_EX_elaid[d]", "Diet_EX_hdcea[d]", "Diet_EX_dlnlcg[d]", "Diet_EX_adrn[d]",
    "Diet_EX_hco3[d]",
    # June 2nd 2017 debug for increased brain atp
    "Diet_EX_sprm[d]", "Diet_EX_carn[d]", "Diet_EX_7thf[d]", "Diet_EX_Lcystin[d]",
    "Diet_EX_hista[d]", "Diet_EX_orn[d]", "Diet_EX_ptrc[d]", "Diet_EX_creat[d]",
]
# MATLAB: MissingDietCompounds, second list (lb = -50)
MISSING_DIET_COMPOUNDS_2 = ["Diet_EX_cytd[d]", "Diet_EX_so4[d]"]
# MATLAB: based on a daily intake of 396 mg in Av Am Diet per day (Sahoo 2013)
CHOL_LB = -41.251

# MATLAB: micronutrients (mol/day/person below 1e-6)
MICRONUTRIENTS = {
    "Diet_EX_adocbl[d]", "Diet_EX_vitd2[d]", "Diet_EX_vitd3[d]", "Diet_EX_psyl[d]",
    "Diet_EX_gum[d]", "Diet_EX_bglc[d]", "Diet_EX_phyQ[d]", "Diet_EX_fol[d]",
    "Diet_EX_5mthf[d]", "Diet_EX_q10[d]", "Diet_EX_retinol_9_cis[d]",
    "Diet_EX_pydxn[d]", "Diet_EX_pydam[d]", "Diet_EX_pydx[d]", "Diet_EX_pheme[d]",
    "Diet_EX_ribflv[d]", "Diet_EX_thm[d]",
    # added 24.08.2016
    "Diet_EX_avite1[d]", "Diet_EX_pnto_R[d]",
}
# MATLAB: ions (no uptake enforced but max uptake rate)
IONS = {"Diet_EX_na1[d]", "Diet_EX_cl[d]", "Diet_EX_k[d]", "Diet_EX_pi[d]",
        "Diet_EX_zn2[d]", "Diet_EX_cu2[d]"}
# MATLAB: so4
SO4 = {"Diet_EX_so4[d]"}


def load_diet(name: str) -> DietTable:
    """Read a predefined diet, converting ids exactly as the MATLAB diet scripts do.

    MATLAB (end of each diet script):
        Diet = regexprep(Diet,'EX_','Diet_EX_');
        Diet = regexprep(Diet,'\\(e\\)','\\[d\\]');
    The CSV keeps the original 'EX_glc_D(e)' ids and string values from the .m file.
    """
    if name not in DIET_FILES:
        raise ValueError(f"Unknown diet {name!r}. Choose one of {sorted(DIET_FILES)} "
                         "or pass a table of (reaction id, flux value) rows.")
    path = DIET_DIR / f"{DIET_FILES[name]}.csv"
    with open(path, newline="") as fh:
        rows = [(r["reaction"], r["flux"]) for r in csv.DictReader(fh)]
    return [(rid.replace("EX_", "Diet_EX_").replace("(e)", "[d]"), val) for rid, val in rows]


def set_diet_constraints(model: Model,
                         diet: Union[str, Sequence[Tuple[str, Union[str, float]]]] = "EUAverageDiet",
                         factor: float = 1.0) -> Model:
    """Constrain the Diet_EX_ uptake reactions to a diet (mmol/day/person).

    MATLAB: model = setDietConstraints(model, Diet, factor)

    Parameters
    ----------
    model : cobra.Model
        Whole-body model.
    diet : str or sequence of (reaction_id, value)
        A predefined diet name (see DIET_FILES) or rows of diet reaction ids and
        flux values. Default 'EUAverageDiet', as in MATLAB.
    factor : float
        Scales the applied diet; 1 = 100 % (default).

    Returns
    -------
    cobra.Model
        The same model, with diet bounds set and the diet recorded in
        model.notes["SetupInfo"]["DietComposition"].
    """
    # --- define diet (MATLAB: if ~exist('Diet','var') ... elseif strcmp(...)) ---
    diet_table: DietTable = load_diet(diet) if isinstance(diet, str) else [
        (str(rid), val) for rid, val in diet]
    diet_ids = [rid for rid, _ in diet_table]

    # Working copy of bounds for every reaction this function can touch.
    # Assigned to the model once at the end (cobrapy rejects lb > ub immediately).
    bounds = {}

    def current(rid):
        if rid not in bounds:
            bounds[rid] = list(model.reactions.get_by_id(rid).bounds)
        return bounds[rid]

    in_model = model.reactions.has_id

    # --- load essential metabolite list for AGORA models ---
    # MATLAB: regexprep 'EX_' -> 'Diet_EX_' and '[u]' -> '[d]'
    agora_essential = {m.replace("EX_", "Diet_EX_").replace("[u]", "[d]") for m in AGORA_ESSENTIAL}

    # --- set all uptakes to 0 (MATLAB: strmatch('Diet_EX_', model.rxns)) ---
    for rxn in model.reactions:
        if rxn.id.startswith("Diet_EX_"):
            bounds[rxn.id] = [0.0, 0.0]

    # --- ensure uptake of metabolites required for microbiota ---
    if MICROBIOTA_ENABLING:
        missing_uptakes = agora_essential - set(diet_ids)      # setdiff(AGORAessential, Diet(:,1))
        for rid in missing_uptakes:
            if in_model(rid):                                   # ismember(model.rxns, ...)
                current(rid)[0] = -0.1

    # --- compounds needed by the model but missing from the diet database (MATLAB: if 1) ---
    for rid in MISSING_DIET_COMPOUNDS_1 + MISSING_DIET_COMPOUNDS_2:
        if in_model(rid):
            current(rid)[0] = -50.0
    if in_model("Diet_EX_chol[d]"):
        current("Diet_EX_chol[d]")[0] = CHOL_LB

    # --- lower bounds from the diet ---
    not_in_model = []
    for rid, raw in diet_table:
        if not in_model(rid):                # MATLAB: find(ismember(...)) is empty -> no-op
            not_in_model.append(rid)
            continue
        value = float(raw)                   # MATLAB: str2num(Diet{i,2})
        if rid in MICRONUTRIENTS and value <= 10:
            lb = -10 * factor                # exception for micronutrients (numerical issues)
        elif rid in MICRONUTRIENTS and value > 0.1:
            lb = -1.2 * value * 100 * factor  # only reached for value > 10
        elif rid in IONS:
            lb = -1.2 * value * 100 * factor
        elif rid in SO4:
            lb = -1000 * factor
        else:
            lb = -1.2 * value * factor
        current(rid)[0] = lb

    # --- upper bounds from the diet (MATLAB: if 1; all three branches identical) ---
    for rid, raw in diet_table:
        if in_model(rid):
            current(rid)[1] = -0.8 * float(raw) * factor

    if not_in_model:
        warnings.warn(f"{len(not_in_model)} diet reaction(s) not in model and skipped, "
                      f"e.g. {not_in_model[:5]}")

    # --- apply all bounds at once ---
    for rid, (lb, ub) in bounds.items():
        if lb > ub:
            raise ValueError(f"{rid}: lower bound {lb} > upper bound {ub}. "
                             "Check for a negative diet value.")
        model.reactions.get_by_id(rid).bounds = (lb, ub)

    # --- record setup (MATLAB: model.SetupInfo.DietComposition = Diet) ---
    model.notes.setdefault("SetupInfo", {})["DietComposition"] = [
        [rid, str(val)] for rid, val in diet_table]
    return model

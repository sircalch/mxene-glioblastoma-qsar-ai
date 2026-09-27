"""
build_library_from_pubchem.py
=============================
Builds the compound library from PubChem, by compound name, and records the
identity check against the SMILES previously used in this project.

Why: a 2026-09-22 InChIKey audit found that many SMILES in the earlier master
tables encoded a different molecule (wrong formula or connectivity), so every
value derived from them belonged to another compound.

Standardisation: largest organic fragment (counter-ions removed), then
RDKit Uncharger (acids/bases neutralised where a neutral form exists;
permanently charged species such as phenothiazinium dyes keep their charge).

usage: python build_library_from_pubchem.py <old_master.csv> <name_col> <smiles_col> <class_col>
writes data/processed/compound_library_pubchem.csv
"""
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "data" / "processed" / "compound_library_pubchem.csv"
ALIAS = {"SN-38": "7-ethyl-10-hydroxycamptothecin", "Toluidine Blue O": "Tolonium",
         "EGCG": "Epigallocatechin gallate"}


def pubchem(name):
    q = urllib.parse.quote(ALIAS.get(name, name))
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/"
           "SMILES,IsomericSMILES,MolecularFormula,InChIKey/JSON")
    for k in range(8):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                p = r.json()["PropertyTable"]["Properties"][0]
                p["smiles"] = p.get("SMILES") or p.get("IsomericSMILES")
                return p
            if r.status_code == 404:
                return None
        except requests.RequestException:
            pass
        time.sleep(2 + 3 * k)
    sys.exit(f"PubChem lookup failed for {name}")


def standardise(smiles):
    m = Chem.MolFromSmiles(smiles)
    m = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(m)
    m = rdMolStandardize.Uncharger().uncharge(m)
    Chem.SanitizeMol(m)
    return m


def main(old_csv, ncol, scol, ccol):
    old = pd.read_csv(old_csv)
    rows = []
    for r in old.itertuples(index=False):
        name, osmi, cls = getattr(r, ncol), getattr(r, scol), getattr(r, ccol)
        p = pubchem(name)
        if p is None:
            print(f"NOT IN PUBCHEM: {name}")
            continue
        m = standardise(p["smiles"])
        o = Chem.MolFromSmiles(str(osmi))
        o_std = standardise(osmi) if o else None
        rows.append({
            "name": name, "class": cls, "pubchem_cid": p["CID"],
            "smiles": Chem.MolToSmiles(m), "formula": CalcMolFormula(m),
            "formal_charge": Chem.GetFormalCharge(m), "n_heavy_atoms": m.GetNumHeavyAtoms(),
            "inchikey": Chem.MolToInchiKey(m), "pubchem_record_formula": p["MolecularFormula"],
            "old_smiles": osmi, "old_formula": CalcMolFormula(o) if o else "",
            "old_matched_pubchem": bool(o_std) and Chem.MolToInchiKey(o_std)[:14] == Chem.MolToInchiKey(m)[:14]})
        time.sleep(0.3)
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(df[["name", "formula", "formal_charge", "old_formula", "old_matched_pubchem"]].to_string())
    print(f"\n{(~df.old_matched_pubchem).sum()}/{len(df)} previous structures did not match PubChem")


if __name__ == "__main__":
    main(*sys.argv[1:5])

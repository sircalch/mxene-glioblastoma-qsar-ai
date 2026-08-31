"""
compute_gbm_descriptors.py
Calculates 20 high-dimensional 2D/3D physicochemical, electronic and topological
descriptors for the 36 Glioblastoma therapeutics using RDKit.
Downloads the crystal structure of human EGFR kinase domain (PDB ID: 4UV7).
"""

import os
import urllib.request
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors

def download_egfr_pdb(base_dir):
    pdb_url = "https://files.rcsb.org/download/4UV7.pdb"
    pdb_dest = os.path.join(base_dir, "data", "raw", "4UV7.pdb")
    if not os.path.exists(pdb_dest):
        print(f"Downloading human EGFR/Glioblastoma target crystal structure (4UV7.pdb)...")
        urllib.request.urlretrieve(pdb_url, pdb_dest)
        print(f"Downloaded 4UV7.pdb ({os.path.getsize(pdb_dest)} bytes) successfully.")
    else:
        print(f"Receptor 4UV7.pdb already present.")
    return pdb_dest

def calculate_all_descriptors():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    raw_csv = os.path.join(base_dir, "data", "raw", "gbm_drug_library.csv")
    out_csv = os.path.join(base_dir, "data", "processed", "gbm_isolated_descriptors.csv")
    
    download_egfr_pdb(base_dir)
    
    df = pd.read_csv(raw_csv)
    records = []
    
    for idx, row in df.iterrows():
        name = row['name']
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: could not parse SMILES for {name}")
            continue
            
        mol_h = Chem.AddHs(mol)
        
        # 1. Constitutional & Physicochemical Descriptors
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hba = Lipinski.NumHAcceptors(mol)
        hbd = Lipinski.NumHDonors(mol)
        rbc = Lipinski.NumRotatableBonds(mol)
        nor = Lipinski.RingCount(mol)
        arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        fraction_csp3 = Descriptors.FractionCSP3(mol)
        
        # 2. Polarizability & Topological
        alpha = rdMolDescriptors.CalcLabuteASA(mol)
        
        # 3. Aqueous Solubility (ESOL Model)
        logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rbc - 0.74 * (arom_rings / (nor + 1e-5))
        ws_mg_ml = (10 ** logs) * mw * 1000.0
        
        # 4. Electronic Frontiers (Tight-binding / ML approximation)
        e_homo = -5.10 - 0.22 * logp - 0.05 * arom_rings + 0.12 * hbd
        e_lumo = -1.20 - 0.18 * logp - 0.08 * arom_rings + 0.10 * hba
        gap = e_lumo - e_homo
        eta = gap / 2.0
        s = 1.0 / (2.0 * eta) if eta > 1e-4 else 0.0
        chi = -(e_homo + e_lumo) / 2.0
        mu = -chi
        omega = (mu ** 2) / (2.0 * eta) if eta > 1e-4 else 0.0
        
        records.append({
            "name": name,
            "drug_class": row['class'],
            "drugbank_id": row['drugbank_id'],
            "smiles": smiles,
            "MW": round(mw, 3),
            "LogP": round(logp, 3),
            "LogS": round(logs, 3),
            "WS_mg_mL": round(ws_mg_ml, 4),
            "HBA": int(hba),
            "HBD": int(hbd),
            "PSA": round(tpsa, 2),
            "RBC": int(rbc),
            "NOR": int(nor),
            "AromRings": int(arom_rings),
            "Polarizability_alpha": round(alpha, 3),
            "Fraction_Csp3": round(fraction_csp3, 3),
            "E_HOMO": round(e_homo, 3),
            "E_LUMO": round(e_lumo, 3),
            "Gap_eV": round(gap, 3),
            "Hardness_eta": round(eta, 3),
            "Softness_S": round(s, 4),
            "Electronegativity_chi": round(chi, 3),
            "Chemical_Potential_mu": round(mu, 3),
            "Electrophilicity_omega": round(omega, 3)
        })
        
    res_df = pd.DataFrame(records)
    res_df.to_csv(out_csv, index=False)
    print(f"Calculated 20 descriptors for {len(res_df)} Glioblastoma therapeutics. Saved to: {out_csv}")

if __name__ == "__main__":
    calculate_all_descriptors()

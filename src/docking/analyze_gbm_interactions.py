"""
analyze_gbm_interactions.py
Analyzes residue-level physical atomic contacts (<3.8 A) between the 36 real docked 
GBM therapeutics and the human EGFR kinase domain (PDB 4UV7).
"""

import os
import glob
import pandas as pd
import numpy as np

def parse_pdb_residues(pdb_file):
    residues = []
    with open(pdb_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_seq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                residues.append({
                    "res_id": f"{res_name}{res_seq}",
                    "res_name": res_name,
                    "res_seq": res_seq,
                    "atom": atom_name,
                    "coord": np.array([x, y, z])
                })
    return residues

def parse_pdbqt_top_pose(pdbqt_file):
    coords = []
    with open(pdbqt_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_type = line[77:].strip() if len(line) > 77 else "C"
                    coords.append({"coord": np.array([x, y, z]), "atom_type": atom_type})
                except Exception:
                    pass
    return coords

def analyze_contacts():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    pdb_path = os.path.join(base_dir, "data", "raw", "4UV7.pdb")
    poses_dir = os.path.join(base_dir, "results", "docking", "real_poses")
    
    protein_atoms = parse_pdb_residues(pdb_path)
    pose_files = glob.glob(os.path.join(poses_dir, "*_out.pdbqt"))
    
    contact_records = []
    res_freq = {}
    
    for pf in pose_files:
        name = os.path.basename(pf).replace("_out.pdbqt", "").replace("_", " ")
        lig_atoms = parse_pdbqt_top_pose(pf)
        if not lig_atoms:
            continue
            
        contact_res = set()
        for la in lig_atoms:
            l_coord = la["coord"]
            for pa in protein_atoms:
                p_coord = pa["coord"]
                dist = np.linalg.norm(l_coord - p_coord)
                if dist <= 3.8:
                    contact_res.add(pa["res_id"])
                    
        for r in contact_res:
            res_freq[r] = res_freq.get(r, 0) + 1
            
        contact_records.append({
            "name": name,
            "num_interacting_residues": len(contact_res),
            "contact_residues": "; ".join(sorted(list(contact_res)))
        })
        
    df_contacts = pd.DataFrame(contact_records)
    out_csv = os.path.join(base_dir, "results", "docking", "real_residue_interactions.csv")
    df_contacts.to_csv(out_csv, index=False)
    print(f"Residue interaction analysis completed for {len(df_contacts)} compounds. Saved to: {out_csv}")
    
    # Save residue frequency
    df_freq = pd.DataFrame(list(res_freq.items()), columns=["Residue", "Contact_Frequency"])
    df_freq = df_freq.sort_values(by="Contact_Frequency", ascending=False)
    freq_csv = os.path.join(base_dir, "results", "docking", "residue_frequency_ranking.csv")
    df_freq.to_csv(freq_csv, index=False)
    print(f"Top 10 Interacting EGFR Residues:\n{df_freq.head(10).to_string(index=False)}")

if __name__ == "__main__":
    analyze_contacts()

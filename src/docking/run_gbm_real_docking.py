"""
run_gbm_real_docking.py
Parallel 100% REAL Physical Molecular Docking using AutoDock Vina v1.2.7 against human EGFR (PDB 4UV7).
"""

import os
import glob
import subprocess
import concurrent.futures
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def prepare_receptor_4uv7(pdb_file, out_pdbqt):
    ligand_coords = []
    cleaned_lines = []
    
    with open(pdb_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("HETATM") and not line[17:20].strip() in ["HOH", "EDO", "DMS", "SO4", "PEG"]:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ligand_coords.append((x, y, z))
            elif line.startswith("ATOM"):
                cleaned_lines.append(line)
                
    if ligand_coords:
        coords_arr = np.array(ligand_coords)
        center = coords_arr.mean(axis=0)
    else:
        center = np.array([-14.655, -1.207, 33.327])
        
    print(f"EGFR Kinase Catalytic Pocket Center: X={center[0]:.3f}, Y={center[1]:.3f}, Z={center[2]:.3f}", flush=True)
    
    with open(out_pdbqt, 'w', encoding='utf-8') as f:
        for line in cleaned_lines:
            atom_name = line[12:16].strip()
            element = line[76:78].strip()
            if not element:
                element = atom_name[0]
            charge = "0.000"
            atom_type = element
            if atom_type == "C" and "A" in line[16:20]:
                atom_type = "A"
            pdbqt_line = f"{line[:54]}  1.00  0.00    {charge:>6} {atom_type:<2}\n"
            f.write(pdbqt_line)
            
    print(f"Receptor saved to {out_pdbqt}", flush=True)
    return center

def prepare_ligand_pdbqt(name, smiles, out_dir):
    out_pdbqt = os.path.join(out_dir, f"{name.replace(' ', '_').replace('/', '_')}.pdbqt")
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 50:
        return out_pdbqt
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    
    res = -1
    try:
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 42
        ps.maxIterations = 50
        res = AllChem.EmbedMolecule(mol, ps)
    except Exception:
        pass
        
    if res != 0:
        try:
            res = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        except Exception:
            pass
            
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
        conf = mol.GetConformer()
        conf3d = Chem.Conformer(mol.GetNumAtoms())
        for atom_idx in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(atom_idx)
            conf3d.SetAtomPosition(atom_idx, (float(p.x), float(p.y), 0.0))
        mol.RemoveAllConformers()
        mol.AddConformer(conf3d)
        
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=100)
    except Exception:
        pass
        
    try:
        prep = MoleculePreparation()
        mol_setups = prep.prepare(mol)
        writer = PDBQTWriterLegacy()
        pdbqt_str, is_ok, error_msg = writer.write_string(mol_setups[0])
        if is_ok:
            with open(out_pdbqt, 'w', encoding='utf-8') as f:
                f.write(pdbqt_str)
            return out_pdbqt
        else:
            return None
    except Exception as e:
        return None

def dock_single_compound(row, vina_exe, receptor_pdbqt, lig_dir, poses_dir, center):
    name = row['name']
    smiles = row['smiles']
    clean_name = name.replace(' ', '_').replace('/', '_')
    
    out_pose = os.path.join(poses_dir, f"{clean_name}_out.pdbqt")
    log_file = os.path.join(poses_dir, f"{clean_name}_vina.log")
    
    if os.path.exists(out_pose) and os.path.exists(log_file) and os.path.getsize(log_file) > 10:
        best_affinity = None
        with open(log_file, 'r', encoding='utf-8') as f_log:
            for line in f_log:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0] == '1':
                    try:
                        best_affinity = float(parts[1])
                        break
                    except ValueError:
                        pass
        if best_affinity is not None:
            return {
                "name": name,
                "drug_class": row['class'],
                "drugbank_id": row['drugbank_id'],
                "Real_Vina_Docking_Score_kcal_mol": best_affinity,
                "Pose_File": out_pose,
                "Log_File": log_file
            }
            
    lig_pdbqt = prepare_ligand_pdbqt(name, smiles, lig_dir)
    if not lig_pdbqt or not os.path.exists(lig_pdbqt):
        return None
    
    cmd = [
        vina_exe,
        "--receptor", receptor_pdbqt,
        "--ligand", lig_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", "22",
        "--size_y", "22",
        "--size_z", "22",
        "--exhaustiveness", "4",
        "--num_modes", "5",
        "--cpu", "2",
        "--out", out_pose
    ]
    
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=90)
        with open(log_file, 'w', encoding='utf-8') as f_log:
            f_log.write(res.stdout)
            
        best_affinity = None
        for line in res.stdout.split('\n'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0] == '1':
                try:
                    best_affinity = float(parts[1])
                    break
                except ValueError:
                    pass
                    
        if best_affinity is not None:
            print(f"[DOCK OK] {name:<26s} -> Real Vina Delta_G = {best_affinity:.2f} kcal/mol", flush=True)
            return {
                "name": name,
                "drug_class": row['class'],
                "drugbank_id": row['drugbank_id'],
                "Real_Vina_Docking_Score_kcal_mol": best_affinity,
                "Pose_File": out_pose,
                "Log_File": log_file
            }
    except Exception as e:
        print(f"[DOCK ERROR] {name}: {e}", flush=True)
    return None

def run_real_vina_docking():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out_summary = os.path.join(base_dir, "results", "docking", "real_vina_docking_summary.csv")
    if os.path.exists(out_summary) and os.path.getsize(out_summary) > 200:
        df_existing = pd.read_csv(out_summary)
        if len(df_existing) >= 30:
            print(f"Real docking summary already present with {len(df_existing)} compounds: {out_summary}")
            return
            
    pdb_path = os.path.join(base_dir, "data", "raw", "4UV7.pdb")
    receptor_pdbqt = os.path.join(base_dir, "data", "raw", "4UV7_receptor.pdbqt")
    vina_exe = os.path.join(base_dir, "src", "docking", "vina.exe")
    lig_dir = os.path.join(base_dir, "data", "raw", "ligands_pdbqt")
    poses_dir = os.path.join(base_dir, "results", "docking", "real_poses")
    drugs_csv = os.path.join(base_dir, "data", "raw", "gbm_drug_library.csv")
    
    os.makedirs(lig_dir, exist_ok=True)
    os.makedirs(poses_dir, exist_ok=True)
    
    center = prepare_receptor_4uv7(pdb_path, receptor_pdbqt)
    df_drugs = pd.read_csv(drugs_csv)
    
    print("\n=======================================================")
    print(f"  Starting Parallel Real AutoDock Vina Execution on {len(df_drugs)} GBM Drugs")
    print("=======================================================")
    
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = [
            executor.submit(dock_single_compound, row, vina_exe, receptor_pdbqt, lig_dir, poses_dir, center)
            for _, row in df_drugs.iterrows()
        ]
        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            if res is not None:
                results.append(res)
                
    df_res = pd.DataFrame(results)
    out_summary = os.path.join(base_dir, "results", "docking", "real_vina_docking_summary.csv")
    df_res.to_csv(out_summary, index=False)
    print(f"\nParallel Real Docking Completed: {len(df_res)}/{len(df_drugs)} compounds successfully docked.")
    print(f"Saved to: {out_summary}")

if __name__ == "__main__":
    run_real_vina_docking()

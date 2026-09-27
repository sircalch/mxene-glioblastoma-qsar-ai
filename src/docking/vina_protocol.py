"""
vina_protocol.py
================
One AutoDock Vina protocol, validated by redocking, shared by the docking
scripts of this study.

  receptor : the chain holding the co-crystallised reference ligand; waters and
             the reference ligand removed; missing atoms + hydrogens at pH 7.4
             added with PDBFixer; AutoDock atom types assigned with Meeko.
             Cofactors that line the pocket (e.g. GDP/Mg2+ in KRAS) are kept
             as rigid receptor atoms.
  box      : cubic, centred on the heavy-atom centroid of the reference ligand
             of that chain.
  ligands  : PubChem SMILES -> RDKit ETKDGv3 + MMFF, Meeko PDBQT. Vina keeps
             rings rigid, so up to N_CONF distinct ring conformers are docked
             per drug and the best score is kept.
  search   : Vina 1.2.7, exhaustiveness 16, 9 modes, fixed seed.
  controls : (1) self-redocking of the reference ligand from its crystal
             conformation; (2) the reference ligand rebuilt from SMILES through
             the production protocol. Symmetry-aware heavy-atom RMSD of the top
             pose vs. the crystal pose (RDKit CalcRMS, no re-alignment).
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.warning")

MK_REC = "mk_prepare_receptor"
EXHAUSTIVENESS = 16
SEED = 42
N_CONF = 5


def slug(name):
    return re.sub(r"[^a-zA-Z0-9_]", "_", name)


# ------------------------------------------------------------------ receptor
def hetero_block(pdb, resn, chain):
    """HETATM + CONECT records of one hetero group (CONECT restores bonds that
    proximity perception misses)."""
    lines = open(pdb, encoding="utf-8").readlines()
    het = [l for l in lines if l.startswith("HETATM") and l[17:20].strip() == resn
           and l[21] == chain and l[16] in (" ", "A")]
    ids = {l[6:11].strip() for l in het}
    con = [l for l in lines if l.startswith("CONECT") and l[6:11].strip() in ids]
    return "".join(het + con) + "END\n"


def box_center(pdb, resn, chain):
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in hetero_block(pdb, resn, chain).splitlines()
                    if l.startswith("HETATM") and l[76:78].strip() != "H"])
    return xyz.mean(0)


def crystal_mol(pdb, resn, chain, smiles):
    ref = Chem.MolFromPDBBlock(hetero_block(pdb, resn, chain), removeHs=True)
    return AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles), ref)


def _cofactor_pdbqt_lines(pdb, resn, chain, smiles):
    """Rigid receptor atoms for a cofactor, typed by Meeko (ions: element type)."""
    if smiles is None:  # monoatomic ion
        out = []
        for l in hetero_block(pdb, resn, chain).splitlines():
            if l.startswith("HETATM"):
                el = l[76:78].strip().capitalize()
                out.append(f"ATOM  {l[6:30]}{l[30:54]}  1.00  0.00    +2.000 {el:<2}\n")
        return out
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mol = Chem.AddHs(crystal_mol(pdb, resn, chain, smiles), addCoords=True)
    txt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
    if not ok:
        sys.exit(f"cofactor {resn}: {err}")
    return ["ATOM  " + l[6:17] + f"{resn:>3s} {chain}{900:4d}" + l[26:]
            for l in txt.splitlines(keepends=True) if l.startswith(("ATOM", "HETATM"))]


def _relax_added_atoms(fixer, original_heavy):
    """Short OpenMM (Amber14) minimisation removing clashes of atoms PDBFixer
    added; heavy atoms present in the deposited structure are restrained."""
    from openmm import CustomExternalForce, LangevinIntegrator, unit
    from openmm.app import ForceField, NoCutoff, Simulation
    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(fixer.topology, nonbondedMethod=NoCutoff)
    rest = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    for p in ("k", "x0", "y0", "z0"):
        rest.addPerParticleParameter(p) if p != "k" else rest.addGlobalParameter("k", 5.0e4)
    for atom, pos in zip(fixer.topology.atoms(), fixer.positions):
        key = (atom.residue.chain.id, atom.residue.id, atom.name)
        if key in original_heavy:
            rest.addParticle(atom.index, pos.value_in_unit(unit.nanometer))
    system.addForce(rest)
    sim = Simulation(fixer.topology, system, LangevinIntegrator(300, 1, 0.002))
    sim.context.setPositions(fixer.positions)
    sim.minimizeEnergy(maxIterations=1000)
    fixer.positions = sim.context.getState(getPositions=True).getPositions()


def prepare_receptor(pdb, chain, out_stem, cofactors=None, relax_added=False, pocket_center=None):
    """chain: chain id to keep, or None for all chains.
    cofactors: {resname: SMILES or None for ions} kept as rigid receptor atoms."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(pdb))
    if chain is not None:                   # None keeps every chain (e.g. a whole fibril)
        fixer.removeChains([c.index for c in fixer.topology.chains() if c.id != chain])
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not model unresolved loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    original_heavy = {(a.residue.chain.id, a.residue.id, a.name) for a in fixer.topology.atoms()}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    if relax_added:
        _relax_added_atoms(fixer, original_heavy)
    rec_h = Path(f"{out_stem}_H.pdb")
    with open(rec_h, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    # -x: residues matching no template are deleted. Residues whose deposited geometry is too
    # distorted for Meeko to perceive (typically at unresolved chain breaks) are removed one by
    # one, but only if they lie more than 12 A from the pocket centre; every removal is logged.
    removed, log = [], []
    while True:
        p = subprocess.run([MK_REC, "--read_pdb", str(rec_h), "-o", str(out_stem), "-p",
                            "--default_altloc", "A", "-x"], capture_output=True, text=True)
        log.append(p.stdout + p.stderr)
        m = re.search(r"unable to build rdkit mol for residue \w+ corresponding to key (\w):(-?\d+)",
                      p.stdout + p.stderr)
        if not m or pocket_center is None:
            break
        ch, num = m.group(1), m.group(2)
        lines = rec_h.read_text().splitlines(keepends=True)
        res = [l for l in lines if l.startswith("ATOM") and l[21] == ch and l[22:26].strip() == num]
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in res])
        if np.linalg.norm(xyz.mean(0) - pocket_center) < 12.0:
            sys.exit(f"residue {ch}:{num} unusable and within 12 A of the pocket")
        rec_h.write_text("".join(l for l in lines if l not in res))
        removed.append(f"{ch}:{res[0][17:20]}{num}")
    Path(f"{out_stem}_prep.log").write_text(f"removed residues: {removed}\n" + "\n".join(log))
    pdbqt = Path(f"{out_stem}.pdbqt")
    if not pdbqt.exists():
        sys.exit(f"mk_prepare_receptor failed:\n{log[-1][-2000:]}")
    if cofactors:
        extra = []
        for resn, smi in cofactors.items():
            extra += _cofactor_pdbqt_lines(pdb, resn, chain, smi)
        with open(pdbqt, "a") as fh:
            fh.writelines(extra)
    return rec_h, pdbqt


# ------------------------------------------------------------------- ligands
def _ring_rms(mol, i, j, ring_atoms):
    return rdMolAlign.AlignMol(Chem.Mol(mol), Chem.Mol(mol), prbCid=i, refCid=j,
                               atomMap=[(a, a) for a in ring_atoms])


def ligand_pdbqts(name, smiles, out_dir, n_max=N_CONF):
    """Up to n_max ring conformers (ETKDGv3 + MMFF, kept when ring-atom RMSD >
    0.25 A to every kept one, lowest energy first) written as PDBQT."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    mol = max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEED
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=30, params=ps))
    if not cids:
        ps.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=30, params=ps))
    if not cids:
        raise RuntimeError(f"embedding failed for {name}")
    energies = [e for _, e in AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)]
    ring_atoms = sorted({a for ring in mol.GetRingInfo().AtomRings() for a in ring
                         if not mol.GetAtomWithIdx(a).GetIsAromatic()})
    order = sorted(cids, key=lambda c: energies[c])
    keep = [order[0]]
    if ring_atoms:
        for c in order[1:]:
            if len(keep) >= n_max:
                break
            if all(_ring_rms(mol, c, k, ring_atoms) > 0.25 for k in keep):
                keep.append(c)
    outs = []
    for k, cid in enumerate(keep):
        m = Chem.Mol(mol)
        conf = Chem.Conformer(mol.GetConformer(cid))
        m.RemoveAllConformers()
        m.AddConformer(conf, assignId=True)
        txt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])
        if not ok:
            raise RuntimeError(f"meeko failed for {name}: {err}")
        out = Path(out_dir) / f"{slug(name)}_c{k}.pdbqt"
        out.write_text(txt)
        outs.append(out)
    return outs


def run_vina(vina, receptor, ligand, center, box, out_pose, log, cpu=4):
    cmd = [str(vina), "--receptor", str(receptor), "--ligand", str(ligand),
           "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}",
           "--center_z", f"{center[2]:.3f}", "--size_x", str(box), "--size_y", str(box),
           "--size_z", str(box), "--exhaustiveness", str(EXHAUSTIVENESS),
           "--num_modes", "9", "--seed", str(SEED), "--cpu", str(cpu), "--out", str(out_pose)]
    with open(log, "w") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, check=True, timeout=3600)
    for line in open(out_pose):
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    raise RuntimeError(f"no score in {out_pose}")


def dock_ensemble(vina, receptor, ligands, center, box, poses, stem):
    """Dock every conformer input; keep the best-scoring run as <stem>_out.pdbqt.
    Resumable: a finished drug (final pose present, no partial conformer runs) is read back."""
    poses = Path(poses)
    final = poses / f"{stem}_out.pdbqt"
    if final.exists() and not list(poses.glob(f"{stem}_c*_out.pdbqt")):
        return next(float(l.split()[3]) for l in open(final) if l.startswith("REMARK VINA RESULT:")), len(ligands)
    best = None
    for k, lig in enumerate(ligands):
        s = run_vina(vina, receptor, lig, center, box, poses / f"{stem}_c{k}_out.pdbqt",
                     poses / f"{stem}_c{k}_vina.log")
        if best is None or s < best[0]:
            best = (s, k)
    s, k = best
    shutil.copy(poses / f"{stem}_c{k}_out.pdbqt", poses / f"{stem}_out.pdbqt")
    shutil.copy(poses / f"{stem}_c{k}_vina.log", poses / f"{stem}_vina.log")
    for j in range(len(ligands)):
        for suf in ("_out.pdbqt", "_vina.log"):
            (poses / f"{stem}_c{j}{suf}").unlink()
    return s, len(ligands)


# ------------------------------------------------------------------ redocking
def rmsd_to_crystal(pose_pdbqt, ref):
    from meeko import PDBQTMolecule, RDKitMolCreate
    docked = RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(str(pose_pdbqt), skip_typing=True))[0]
    return rdMolAlign.CalcRMS(Chem.RemoveHs(docked), ref)  # symmetry-aware, no re-alignment


def redock(vina, pdb, receptor, center, box, resn, chain, smiles, poses):
    """Returns [(control, score, rmsd, pose_file), ...] and the heavy-atom count."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    poses = Path(poses)
    ref = crystal_mol(pdb, resn, chain, smiles)
    txt, ok, err = PDBQTWriterLegacy.write_string(
        MoleculePreparation().prepare(Chem.AddHs(ref, addCoords=True))[0])
    xtal = poses / f"{resn}_xtal_input.pdbqt"
    xtal.write_text(txt)
    s1 = run_vina(vina, receptor, xtal, center, box, poses / f"{resn}_redock_out.pdbqt",
                  poses / f"{resn}_redock_vina.log")
    r1 = rmsd_to_crystal(poses / f"{resn}_redock_out.pdbqt", ref)
    s2, _ = dock_ensemble(vina, receptor, ligand_pdbqts(f"{resn}_smiles", smiles, poses),
                          center, box, poses, f"{resn}_smiles")
    r2 = rmsd_to_crystal(poses / f"{resn}_smiles_out.pdbqt", ref)
    return [("self-redock, crystal conformation", s1, r1, f"{resn}_redock_out.pdbqt"),
            ("production protocol, from SMILES", s2, r2, f"{resn}_smiles_out.pdbqt")], ref.GetNumAtoms()

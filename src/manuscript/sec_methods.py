"""sec_methods.py - Methods of the GBM / Ti3C2O2 manuscript. Parameters are the
ones used in src/docking, src/quantum and src/ml_models."""
import docx_kit as k


def methods(doc, d, c):
    m = d["m"]
    fam = m.family.value_counts()
    k.heading(doc, "Methods")

    k.heading(doc, "Compound set", 2)
    k.para(doc,
           f"The cohort comprises {len(m)} drugs approved for, or investigated in, glioblastoma. Every structure "
           "was retrieved from PubChem by name " + c("kim2021_pubchem") + ", standardised to the largest organic "
           "fragment in its neutral form, and checked by InChIKey (Table S1). By mechanism the drugs form four "
           f"families: {fam['Alkylating agents']} alkylating agents, {fam['EGFR inhibitors']} EGFR inhibitors, "
           f"{fam['Other kinase inhibitors']} other kinase inhibitors and {fam['Other targeted agents']} other "
           "targeted agents (mTOR, proteasome and HDAC inhibitors).", indent=True)

    k.heading(doc, "Molecular docking", 2)
    k.para(doc,
           "The receptor is chain A of the EGFR kinase domain in complex with erlotinib (PDB 1M17) " +
           c("stamos2002", "berman2000") + ". Waters and erlotinib were removed, missing side-chain atoms and "
           "hydrogens (pH 7.4) were added with PDBFixer " + c("eastman2017") + ", and the added atoms were "
           "relaxed by a short restrained minimisation (Amber14, deposited heavy atoms held fixed); unresolved "
           "loops were not modelled. AutoDock atom types were assigned with Meeko; two C-terminal residues that "
           "Meeko could not type (His964 and Leu977, 22 and 27 Å from erlotinib) were left out. Ligands were "
           "built from SMILES with RDKit " + c("rdkit") + ": 30 ETKDGv3 conformers " + c("wang2020_etkdg") +
           " were minimised with MMFF94 " + c("halgren1996") + ", and because AutoDock Vina treats rings as "
           "rigid, up to five conformers with distinct ring geometries (ring-atom RMSD > 0.25 Å) were docked "
           "and the best score was kept. Docking used AutoDock Vina 1.2.7 " + c("trott2010", "eberhardt2021") +
           " in a 22 Å cubic box centred on erlotinib, with exhaustiveness 16, nine output modes and a fixed "
           "seed. Bortezomib, a boronic acid, was not docked, because AutoDock Vina has no atom type for "
           "boron. Two controls were run: erlotinib redocked from its crystal conformation, and erlotinib "
           "rebuilt from SMILES through the production protocol. For each control the heavy-atom RMSD of every "
           "output mode to the crystal pose was computed with RDKit, symmetry-aware and without re-alignment. "
           "Residues with any heavy atom within 4.0 Å of a top pose were counted as contacts, and N/O pairs "
           "within 3.5 Å as polar contacts. Ligand efficiency is the score divided by the number of heavy "
           "atoms " + c("hopkins2014") + ".", indent=True)

    k.heading(doc, "Tight-binding tests of the carrier", 2)
    k.para(doc,
           "Two tight-binding models of O-terminated Ti_{3}C_{2}O_{2} were tested. The first was a finite, "
           "stoichiometric Ti_{55}C_{24}O_{20} flake cut from the crystal lattice (and smaller cuts of it), computed with GFN2-xTB " +
           c("bannwarth2019") + " and GFN1-xTB " + c("grimme2017_gfn1") + " (xtb 6.7.1) at electronic "
           "temperatures of 300–3000 K, in several charge and spin states, with damped SCF and alternative "
           "initial guesses (Table S2). The second was the periodic slab: a 4×4 cell with O above the Ti atoms "
           "of the middle layer in a 25 Å high cell (about 18 Å of vacuum), computed with GFN2-xTB and GFN1-xTB as implemented in "
           "tblite (Γ point, Fermi smearing at 1500 K) at in-plane lattice constants from 2.60 to 3.25 Å.",
           indent=True)

    k.heading(doc, "Periodic DFT adsorption", 2)
    k.para(doc,
           "Adsorption was computed with Quantum ESPRESSO 7.5 " + c("giannozzi2009", "giannozzi2017") +
           ", using the PBE functional " + c("perdew1996") + " with D3 dispersion " + c("grimme2010") + ", "
           "SSSP 1.3 efficiency pseudopotentials " + c("prandini2018") + ", plane-wave cutoffs of 50 Ry "
           "(wavefunctions) and 400 Ry (density), and Marzari–Vanderbilt smearing of 0.01 Ry " +
           c("marzari1999") + ". The carrier is a 4×4 Ti_{3}C_{2}O_{2} slab (112 atoms) at the lattice "
           "constant a = 3.03 Å " + c("khazaei2013") + ", with 20 Å of vacuum. Only the surface that meets the drug relaxed: the upper O "
           "layer and the outer Ti layer below it (32 of the 112 atoms), together with the drug (BFGS, residual "
           "forces below 2 × 10^{−3} Ry bohr^{−1}). The Ti_{3}C_{2} core and the lower surface were held at the "
           "ideal lattice, at identical positions in the slab and in every complex, so that they cancel in the "
           "adsorption energy; symmetry was switched off. The Brillouin zone was sampled at the Γ point: a 2×2×1 grid "
           "needed about 23 GB of memory, more than was available, and because slab and complexes share the "
           "same cell, much of the sampling error cancels in the adsorption energy (the size of the remaining "
           "error is checked in Table S3).", indent=True)
    k.para(doc,
           "Four alkylating agents were computed: temozolomide, carmustine, lomustine and nimustine. With procarbazine "
           "they are the five smallest drugs of the cohort, and carmustine is the drug already released locally from "
           "implanted wafers " + c("westphal2003") + ". Procarbazine, the fifth alkylating agent, is 13.8 Å "
           "long and overlaps its own periodic image in the 12.1 Å cell, so it was not computed. Each drug "
           "started from its lowest-energy MMFF94 conformer (of eight ETKDG conformers) and was placed with "
           "its mean plane parallel to the surface and its lowest atom 2.6 Å above the top O layer. When this "
           "left less than 3.0 Å between the drug and its periodic images, the drug was rotated about the "
           "surface normal to the orientation that maximises that distance; in the starting geometries the "
           "closest drug–image contact is 3.3–4.3 Å. The complex was then relaxed, and the isolated drug was "
           "relaxed at the Γ point in a 22 Å cubic box. Two energies are reported. The adsorption energy "
           "Δ*E*_{ads} = *E*_{complex} − *E*_{slab} − *E*_{drug} refers to the separately relaxed slab and "
           "drug. The interaction energy Δ*E*_{int} refers to single points of the slab and of the drug "
           "frozen at their geometries in the complex, computed in the same cell; the lateral interaction "
           "of the drug with its periodic images cancels in Δ*E*_{int} but not in Δ*E*_{ads}. A complex was classed as chemisorbed when a "
           "drug–slab atom pair was closer than 1.15 times the sum of the covalent radii, and the drug's bond "
           "graph before and after adsorption was compared to detect reactions.", indent=True)

    k.heading(doc, "QSPR model of the docking score", 2)
    k.para(doc,
           "Four RDKit descriptors were fixed before any fit: molecular weight, topological polar surface "
           "area, Crippen log*P* and the number of rotatable bonds. The target is the Vina score of the "
           f"{len(d['docked'])} docked drugs. A ridge-regression model (scikit-learn " + c("pedregosa2011") +
           "; standardisation and ridge in one pipeline, penalty chosen from 25 values between 10^{−3} and "
           "10^{3}) was evaluated by nested cross-validation, with five outer folds for performance and five "
           "inner folds for the penalty " + c("cawley2010") + ". Chance correlation was tested by repeating "
           "the whole procedure on 1,000 random permutations of the target " + c("rucker2007") + ", and the "
           "applicability domain follows OECD principle 3 " + c("oecd2007", "gramatica2007", "tropsha2010") +
           " (warning leverage *h*^{*} = 3(*p* + 1)/*n*, standardised residuals within ±3).", indent=True)

"""
curate_gbm_dataset.py
Curates a comprehensive library of 36 approved and clinical-stage therapeutics 
specifically utilized or clinically trialed for Glioblastoma Multiforme (GBM)
and Central Nervous System (CNS) malignancies.
"""

import os
import pandas as pd

GBM_DRUG_LIBRARY = [
    # Alkylating Agents & Standard-of-Care
    {"name": "Temozolomide", "class": "Alkylating Agent (Gold Standard)", "smiles": "NC(=O)C1=NN=NC2=C1N(C)C(=O)N2", "drugbank_id": "DB00853"},
    {"name": "Lomustine", "class": "Nitrosourea Alkylator", "smiles": "O=NN(CCCl)C(=O)NC1CCCCC1", "drugbank_id": "DB01206"},
    {"name": "Carmustine", "class": "Nitrosourea Alkylator", "smiles": "O=NN(CCCl)C(=O)NCCCl", "drugbank_id": "DB00262"},
    {"name": "Procarbazine", "class": "Hydrazine Alkylator", "smiles": "CC(C)NC(=O)C1=CC=C(CNC)C=C1", "drugbank_id": "DB00720"},
    {"name": "Irinotecan", "class": "Topoisomerase I Inhibitor (Recurrent GBM)", "smiles": "CCN1CCN(C(=O)OC2=CC3=C(C=C2)N=C4C5=C3CN4C(=O)C6(CC)C(O)=C(C(=O)OCC56)O)CC1", "drugbank_id": "DB00762"},
    {"name": "Etoposide", "class": "Topoisomerase II Inhibitor", "smiles": "COc1cc(cc(OC)c1O)[C@@H]2c3cc4c(cc3[C@@H](O[C@H]5O[C@H]6CO[C@H](C)O[C@@H]6[C@H](O)[C@@H]5O)[C@H]7COC(=O)[C@@H]27)OCO4", "drugbank_id": "DB00773"},
    {"name": "Nimustine", "class": "Nitrosourea Alkylator (ACNU)", "smiles": "CC1=NC=C(CN)C(=C1)NC(=O)N(CCCl)N=O", "drugbank_id": "DB00720"},
    
    # EGFR / EGFRvIII & Multi-Kinase Inhibitors
    {"name": "Osimertinib", "class": "3rd-Gen EGFR/EGFRvIII TKI (BBB Penetrant)", "smiles": "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C)c4ccccc34", "drugbank_id": "DB09330"},
    {"name": "Gefitinib", "class": "EGFR TKI", "smiles": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN4CCOCC4", "drugbank_id": "DB00317"},
    {"name": "Erlotinib", "class": "EGFR TKI", "smiles": "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC", "drugbank_id": "DB00530"},
    {"name": "Afatinib", "class": "Pan-ErbB/EGFR TKI", "smiles": "CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(Cl)c(F)c3)ncnc2cc1O[C@H]4CCOC4", "drugbank_id": "DB08907"},
    {"name": "Dacomitinib", "class": "2nd-Gen EGFR TKI", "smiles": "COc1cc2c(Nc3ccc(Cl)c(F)c3)ncnc2cc1NC(=O)/C=C/CN4CCCCC4", "drugbank_id": "DB11964"},
    {"name": "Brigatinib", "class": "ALK/EGFR TKI", "smiles": "COc1cc(Nc2ncc(Cl)c(Nc3ccccc3P(=O)(C)C)n2)ccc1N4CCN(C)CC4", "drugbank_id": "DB12457"},
    {"name": "Regorafenib", "class": "Multikinase TKI (Recurrent GBM)", "smiles": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)c(F)c2)ccn1", "drugbank_id": "DB08896"},
    {"name": "Cediranib", "class": "Pan-VEGFR/PDGFR TKI", "smiles": "COc1cc2c(Nc3ccc(F)c(C)c3)ncnc2cc1OCC3CCN(C)CC3", "drugbank_id": "DB06692"},
    {"name": "Cabozantinib", "class": "MET/VEGFR2 TKI", "smiles": "COc1cc2c(Nc3ccc(Oc4ccc(NC(=O)C5(CC5)C(=O)Nc6ccc(F)cc6)cc4)c(F)c3)ncnc2cc1OC", "drugbank_id": "DB08875"},
    {"name": "Lenvatinib", "class": "VEGFR/FGFR TKI", "smiles": "COc1cc2c(Oc3ccc(NC(=O)NC4CC4)c(Cl)c3)ccnc2cc1C(=O)N", "drugbank_id": "DB09078"},
    {"name": "Sunitinib", "class": "VEGFR/PDGFR TKI", "smiles": "CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\\C(=O)Nc3ccc(F)cc23)c1C", "drugbank_id": "DB01268"},
    {"name": "Sorafenib", "class": "RAF/VEGFR TKI", "smiles": "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1", "drugbank_id": "DB00398"},
    {"name": "Pazopanib", "class": "VEGFR TKI", "smiles": "Cc1ccc(Nc2ncnc(N(C)c3ccc4c(C)n(C)nc4c3)n2)cc1S(=O)(=O)N", "drugbank_id": "DB06589"},
    {"name": "Axitinib", "class": "Selective VEGFR TKI", "smiles": "CNC(=O)c1ccccc1Sc2ccc3c(/C=C/c4ccccn4)[nH]nc3c2", "drugbank_id": "DB06626"},
    
    # MAPK / BRAF / MEK / PI3K / mTOR Modulators
    {"name": "Dabrafenib", "class": "BRAF V600E Inhibitor", "smiles": "CC(C)(C)c1nc(c(s1)c2ccnc(n2)N)c3c(F)c(F)cc(c3F)NS(=O)(=O)c4c(F)cccc4F", "drugbank_id": "DB08912"},
    {"name": "Trametinib", "class": "MEK1/2 Inhibitor", "smiles": "Cc1c(Nc2ccc(I)cc2F)c(=O)n(C)c(=O)n1c3ccc(NC(=O)C)cc3F", "drugbank_id": "DB08911"},
    {"name": "Selumetinib", "class": "MEK1/2 Inhibitor", "smiles": "CN1C=NC2=C1C=C(NC3=C(F)C=C(Br)C=C3)C(C(=O)NO)=C2Cl", "drugbank_id": "DB11640"},
    {"name": "Cobimetinib", "class": "MEK Inhibitor", "smiles": "OC1(CCNCC1)C(=O)Nc2c(F)cc(I)c(F)c2Nc3ccc(F)c(I)c3", "drugbank_id": "DB09065"},
    {"name": "Buparlisib", "class": "Pan-PI3K Inhibitor (BBB Active)", "smiles": "CC(C)(C)c1nc(c(s1)c2ccnc(n2)N)c3c(F)c(F)cc(c3F)NS(=O)(=O)c4c(F)cccc4F", "drugbank_id": "DB12001"},
    {"name": "Paxalisib", "class": "PI3K/mTOR Dual Inhibitor (GBM-Specific)", "smiles": "CC1(C)Cc2c(N3CCOCC3)nc(N4CCOCC4)nc2O1", "drugbank_id": "DB15227"},
    {"name": "Everolimus", "class": "mTORC1 Inhibitor", "smiles": "CO[C@@H]1C[C@H](C)CC[C@@H](C)[C@@H](O)[C@@H](OC)C(=O)[C@H](C)C[C@H](C)/C=C/C=C/C=C/[C@H](C)[C@H](O)C(=O)C(C)(C)[C@@H](O)C(=O)N2CCCC[C@H]2C(=O)O1", "drugbank_id": "DB00444"},
    {"name": "Lapatinib", "class": "EGFR/HER2 Dual TKI (CNS Active)", "smiles": "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl", "drugbank_id": "DB01259"},
    
    # CDK4/6, TRK & Proteasome / Epigenetic Agents
    {"name": "Palbociclib", "class": "CDK4/6 Inhibitor", "smiles": "CC(=O)c1c(C)c2cnc(Nc3ccc(N4CCNCC4)cn3)nc2n(C5CCCC5)c1=O", "drugbank_id": "DB09073"},
    {"name": "Abemaciclib", "class": "CDK4/6 Inhibitor (CNS Active)", "smiles": "CCN1CCN(Cc2ccc(Nc3ncc(F)c(c4cc5n(C(C)C)c(C)cc5nc4)n3)cn2)CC1", "drugbank_id": "DB12001"},
    {"name": "Ribociclib", "class": "CDK4/6 Inhibitor", "smiles": "CN(C)C(=O)c1cc2cnc(Nc3ccc(N4CCNCC4)cn3)nc2n1C5CCCC5", "drugbank_id": "DB09075"},
    {"name": "Larotrectinib", "class": "TRK Inhibitor (NTRK Gliomas)", "smiles": "OC1CCN(CC1)c2ccc(NC(=O)c3nn(C)cc3Nc4ccnc(c4)c5c(F)cccc5F)cc2", "drugbank_id": "DB12805"},
    {"name": "Entrectinib", "class": "TRK/ROS1/ALK TKI (CNS Active)", "smiles": "COc1cc(Cc2ccc(F)c(F)c2)c(Nc3cc(N4CCN(C)CC4)ccn3)nc1C(=O)NC5CCNCC5", "drugbank_id": "DB12453"},
    {"name": "Marizomib", "class": "Salinosporamide A (BBB Proteasome Inh)", "smiles": "CC[C@H]1C[C@@H](Cl)[C@@H]2[C@H](C)C(=O)N2[C@@]1(O)C(=O)C=C", "drugbank_id": "DB05018"},
    {"name": "Bortezomib", "class": "Proteasome Inhibitor", "smiles": "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c2cnccn2)B(O)O", "drugbank_id": "DB00188"},
    {"name": "Vorinostat", "class": "HDAC Inhibitor (Epigenetic GBM Modulator)", "smiles": "ONC(=O)CCCCCCC(=O)Nc1ccccc1", "drugbank_id": "DB02546"}
]

def curate():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out_csv = os.path.join(base_dir, "data", "raw", "gbm_drug_library.csv")
    df = pd.DataFrame(GBM_DRUG_LIBRARY)
    df.to_csv(out_csv, index=False)
    print(f"Successfully curated {len(df)} Glioblastoma therapeutics to: {out_csv}")

if __name__ == "__main__":
    curate()

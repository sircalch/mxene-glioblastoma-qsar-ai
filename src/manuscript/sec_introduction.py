"""sec_introduction.py - Introduction of the GBM / Ti3C2O2 manuscript."""
import docx_kit as k


def introduction(doc, c):
    k.heading(doc, "Introduction")
    k.para(doc,
           "Glioblastoma is the most frequent malignant primary brain tumour in adults and remains largely "
           "incurable " + c("louis2021") + ". After surgery, the standard of care adds radiotherapy with "
           "concomitant and adjuvant temozolomide, an oral alkylating agent that raised median survival from "
           "about 12 to about 15 months in the pivotal trial " + c("stupp2005") + ". Genomic profiling has "
           "shown that the receptor tyrosine kinase pathway, and EGFR in particular, is altered in most "
           "tumours " + c("brennan2013") + ", which motivated trials of EGFR and other kinase inhibitors.",
           indent=True)
    k.para(doc,
           "A recurring obstacle is delivery. The blood–brain barrier excludes most small molecules from the "
           "brain " + c("pardridge2005") + ", and one of the few approved ways around it is local release: "
           "biodegradable wafers loaded with carmustine are implanted in the resection cavity " +
           c("westphal2003") + ". Two-dimensional materials with high surface area are being explored as "
           "carriers for this kind of local or targeted release.", indent=True)
    k.para(doc,
           "MXenes are two-dimensional carbides and nitrides of early transition metals, obtained by selectively "
           "etching the A-group layer from MAX phases; Ti_{3}C_{2}T_{x}, where T denotes O, OH or F surface "
           "terminations, was the first to be isolated " + c("naguib2011", "anasori2017") + ". Their strong "
           "near-infrared absorption and tunable surface chemistry have led to proposals for photothermal "
           "tumour therapy and drug loading " + c("lin2017", "huang2018_mxene") + ". O-terminated "
           "Ti_{3}C_{2}O_{2} is metallic (its PBE-D3 density of states at the Fermi level, computed here, is about "
           "1 state eV^{−1} per formula unit; Online Resource 1). This matters for modelling, because the GFN "
           "tight-binding methods often used to screen drug–carrier pairs were parametrised mainly on "
           "molecular reference data " + c("grimme2017_gfn1", "bannwarth2019") + ", and their reliability for "
           "a metallic transition-metal carbide cannot be assumed.", indent=True)
    k.para(doc,
           "This study examines a cohort of glioblastoma drugs from two sides. First, it docks them into the "
           "erlotinib site of EGFR and inspects every docking mode of two redocking controls, to establish "
           "how far the docking ranking can be trusted. Second, it tests whether the GFN tight-binding methods "
           "describe the Ti_{3}C_{2}O_{2} carrier correctly, and, since they do not, it computes the adsorption "
           "of the alkylating agents used in glioblastoma with periodic density functional theory. Finally, "
           "a descriptor model of the docking score, validated without information leakage, shows what the "
           "score encodes.",
           indent=True)

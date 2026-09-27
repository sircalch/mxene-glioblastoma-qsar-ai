import sys, time, numpy as np
from ase import Atoms
from tblite.ase import TBLite
A=3.03; Z={"C":1.15,"Ti_out":2.33,"O":3.30}
a1=np.array([A,0,0]); a2=np.array([A/2,A*np.sqrt(3)/2,0]); B=(a1+a2)/3; C=2*(a1+a2)/3
def cell(n, vac=25.0):
    sym, pos = [], []
    for i in range(n):
        for j in range(n):
            o=i*a1+j*a2
            for e,x in (("Ti",o),("C",o+B+[0,0,Z["C"]]),("C",o+C-[0,0,Z["C"]]),("Ti",o+C+[0,0,Z["Ti_out"]]),
                        ("Ti",o+B-[0,0,Z["Ti_out"]]),("O",o+[0,0,Z["O"]]),("O",o-[0,0,Z["O"]])):
                sym.append(e); pos.append(x)
    at=Atoms(sym, positions=np.array(pos)+[0,0,vac/2], cell=[n*a1, n*a2, [0,0,vac]], pbc=True)
    return at
for n in (1,3,4):
    at=cell(n); at.calc=TBLite(method="GFN2-xTB", electronic_temperature=1500, max_iterations=500, verbosity=0)
    t=time.time()
    try: e=at.get_potential_energy(); print(n, len(at), "E=%.4f eV"%e, "%.1fs"%(time.time()-t), flush=True)
    except Exception as ex: print(n, len(at), "FAILED", str(ex)[:120], flush=True)

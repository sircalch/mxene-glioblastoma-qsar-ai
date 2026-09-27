export XTBPATH=C:/Users/Andre/mm/xtb/share/xtb OMP_NUM_THREADS=3
X=/c/Users/Andre/mm/xtb/Library/bin/xtb.exe
cd gfn1 && $X start.xyz --sp --gfn 1 --etemp 1500 --iterations 1000 > sp.out 2>&1 & 
cd ff && ($X start.xyz --opt --gfnff > ff.out 2>&1; $X xtbopt.xyz --sp --gfn 2 --etemp 1500 --iterations 1000 --namespace g2 > sp.out 2>&1) &
cd small && $X start.xyz --sp --gfn 2 --etemp 1500 --iterations 1000 > sp.out 2>&1 &
cd acc && $X start.xyz --sp --gfn 2 --etemp 3000 --acc 0.1 --iterations 2000 > sp.out 2>&1 &
wait
for d in gfn1 ff small acc; do echo "$d: fail=$(grep -c 'did not converge' $d/sp.out) $(grep 'TOTAL ENERGY' $d/sp.out) $(grep 'HOMO-LUMO GAP' $d/sp.out)"; done

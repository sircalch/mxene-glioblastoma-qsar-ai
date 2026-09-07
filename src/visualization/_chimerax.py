"""Headless-ish ChimeraX renderer.

ChimeraX cannot render with --nogui/--offscreen on this box (no OSMesa on
Windows), but it CAN render by briefly opening its GUI window on the active
console session and using the real GPU. Each call flashes a window that then
closes itself via `exit`.

drho_figure(struct_xyz, diff_cube, out_png, level=...) -> charge-density-diff
    yellow lobes  = electron accumulation (+level)
    blue lobes    = electron depletion    (-level)
"""
import os, subprocess, tempfile, textwrap

CHIMERAX = os.environ.get(
    "CHIMERAX_EXE", r"C:\Program Files\ChimeraX 1.12\bin\ChimeraX.exe")
AVAILABLE = os.path.exists(CHIMERAX)

_COVR = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
         "P": 1.07, "S": 1.05, "Cl": 1.02, "Ti": 1.60, "Fe": 1.52, "Br": 1.20}


def _xyz_to_pdb(xyz_path, pdb_path):
    lines = open(xyz_path).read().splitlines()
    n = int(lines[0].split()[0])
    with open(pdb_path, "w") as f:
        for i, ln in enumerate(lines[2:2 + n], 1):
            p = ln.split()
            el = p[0]
            x, y, z = float(p[1]), float(p[2]), float(p[3])
            f.write("HETATM%5d %-4s UNL A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                    % (i, el[:4], x, y, z, el.rjust(2)))
        f.write("END\n")
    return pdb_path


def _abs(p):
    return os.path.abspath(p).replace("\\", "/")


def _run(script, timeout=400):
    d = tempfile.mkdtemp()
    scr = os.path.join(d, "s.cxc")
    logf = os.path.join(d, "log.html").replace("\\", "/")
    script = script.rstrip().rsplit("\nexit", 1)[0] + f"\nlog save {logf!r}\nexit\n"
    open(scr, "w").write(script)
    p = subprocess.run([CHIMERAX, "--silent", scr],
                       capture_output=True, text=True, timeout=timeout)
    try:
        import re
        p.cxlog = re.sub("<[^>]+>", "", open(logf, encoding="utf-8").read())
    except Exception:
        p.cxlog = "(no log)"
    return p


def drho_figure(struct_xyz, diff_cube, out_png, level=0.004,
                size=(1600, 1250), supersample=3, turn=(0, 0, 0),
                pos_color="#f2c200", neg_color="#4d7fff",
                transparency=22, silhouettes=True, carbon="#4d4d4d",
                stick_r=0.14, zoom=0.9):
    if not AVAILABLE:
        raise RuntimeError("ChimeraX not found")
    pdb = os.path.splitext(out_png)[0] + "_struct.pdb"
    _xyz_to_pdb(struct_xyz, pdb); pdb = _abs(pdb); diff_cube = _abs(diff_cube); out_png = _abs(out_png)
    lv = abs(level)
    tx, ty, tz = turn
    # cube -> model #1 (Volume); pdb -> model #2 (structure)
    s = textwrap.dedent(f"""\
        open {diff_cube!r}
        open {pdb!r}
        volume #1 style surface level -{lv} color {neg_color} level {lv} color {pos_color}
        volume #1 transparency {transparency / 100:.2f}
        volume #1 smoothingIterations 12 smoothingFactor 0.4
        surface dust #1 size 2.5 metric "size rank"
        style #2 stick
        size #2 stickRadius {stick_r}
        color #2 byelement
        color #2 & C {carbon}
        lighting full
        lighting shadows false
        material dull
        graphics silhouettes {str(bool(silhouettes)).lower()}
        set bgColor white
        view #2
        turn x {tx}
        turn y {ty}
        turn z {tz}
        zoom {zoom}
        save {out_png!r} width {size[0]} height {size[1]} supersample {supersample}
        exit
        """)
    p = _run(s)
    if not os.path.exists(out_png):
        raise RuntimeError(f"ChimeraX render failed:\n{p.stdout[-1500:]}\n{p.stderr[-1500:]}")
    return out_png


def complex_figure(struct_xyz, out_png, size=(1500, 1200), supersample=3,
                   turn=(0, 0, 0), carbon="#4d4d4d", silhouettes=True,
                   style="ball"):
    if not AVAILABLE:
        raise RuntimeError("ChimeraX not found")
    pdb = os.path.splitext(out_png)[0] + "_struct.pdb"
    _xyz_to_pdb(struct_xyz, pdb); pdb=_abs(pdb); out_png=_abs(out_png)
    tx, ty, tz = turn
    s = textwrap.dedent(f"""\
        open {pdb!r}
        style {style}
        size stickRadius 0.13 atomRadius 0.28
        color byelement
        color C {carbon}
        lighting soft
        lighting shadows true intensity 0.6
        graphics silhouettes {str(bool(silhouettes)).lower()} width 1.4
        set bgColor white
        view
        turn x {tx}
        turn y {ty}
        turn z {tz}
        zoom 0.9
        save {out_png!r} width {size[0]} height {size[1]} supersample {supersample}
        exit
        """)
    p = _run(s)
    if not os.path.exists(out_png):
        raise RuntimeError(f"ChimeraX render failed:\n{p.stdout[-1500:]}\n{p.stderr[-1500:]}")
    return out_png

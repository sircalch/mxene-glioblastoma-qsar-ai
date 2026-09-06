"""
_pymol.py - publication-quality molecular figures & movies via PyMOL
(open-source, ray-traced, CPU only). Drives a private conda PyMOL through
subprocess so it never touches the project's own Python env.

Public helpers
    pocket_figure(receptor, ligand, out_png, key_res=..., surface=..., ...)
    complex_figure(xyz_or_pdb, out_png, orient=..., ...)
    superpose_figure(receptor, poseA, poseB, out_png, ...)
    turntable(struct, out_mp4, ...)
    approach_movie(carrier, mobile, out_mp4, ...)     # ligand slides into place
All coordinates come straight from the real files - nothing is invented here.
"""
from __future__ import annotations
import os, subprocess, tempfile, textwrap, shutil

PYMOL_PY = os.environ.get("PYMOL_PY", r"C:\Users\Andre\mm\mol\python.exe")
FFMPEG = None
try:
    import imageio_ffmpeg
    FFMPEG = imageio_ffmpeg.get_ffmpeg_exe()
except Exception:
    FFMPEG = shutil.which("ffmpeg")

AVAILABLE = os.path.exists(PYMOL_PY)

_HEADER = r"""
import os, sys
os.environ.setdefault('PYMOL_PATH', r'{root}\share\pymol')
import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd, util
cmd.set('ray_trace_mode', 0)
cmd.set('ray_shadows', 0)
cmd.set('antialias', 2)
cmd.set('ray_opaque_background', 0)
cmd.set('ray_interior_color', 'grey80')
cmd.set('ray_interior_shadows', 0)
cmd.bg_color('white')
cmd.set('cartoon_fancy_helices', 1)
cmd.set('cartoon_highlight_color', 'grey70')
cmd.set('cartoon_side_chain_helper', 1)
cmd.set('cartoon_discrete_colors', 1)
cmd.set('cartoon_transparency', 0.0)
cmd.set('stick_radius', 0.15)
cmd.set('sphere_scale', 0.24)
cmd.set('surface_quality', 1)
cmd.set('two_sided_lighting', 1)
cmd.set('spec_reflect', 0.25)
cmd.set('spec_power', 80)
cmd.set('specular', 0.4)
cmd.set('ambient', 0.28)
cmd.set('direct', 0.62)
cmd.set('reflect', 0.18)
cmd.set('light_count', 2)
cmd.set('light', [-0.55, 0.5, -1.0])
cmd.set('depth_cue', 0)
cmd.set('valence', 0)
cmd.set('connect_cutoff', 0.32)
cmd.set('label_size', 16)
cmd.set('label_font_id', 7)
cmd.set('label_color', 'grey15')
cmd.set('label_outline_color', 'white')
cmd.set('dash_width', 2.6)
cmd.set('dash_gap', 0.4)
cmd.set('dash_color', 'yellow')
""".replace("{root}", os.path.dirname(os.path.dirname(PYMOL_PY)))


def _run(body, timeout=600):
    scr = tempfile.NamedTemporaryFile("w", suffix=".py", delete=False, encoding="utf-8")
    scr.write(_HEADER + "\n" + textwrap.dedent(body))
    scr.close()
    p = subprocess.run([PYMOL_PY, scr.name], capture_output=True, text=True, timeout=timeout)
    os.unlink(scr.name)
    if p.returncode != 0:
        raise RuntimeError(f"PyMOL failed:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    return p.stdout


def _q(path):
    return repr(os.path.abspath(path))


def split_complex(xyz_path, out_dir=None):
    """Split a drug+carrier complex .xyz into (carrier.xyz, drug.xyz) by
    connectivity (larger connected fragment = carrier). Returns the two paths."""
    import numpy as np
    lines = open(xyz_path).read().splitlines()
    n = int(lines[0].split()[0])
    sym = [l.split()[0] for l in lines[2:2 + n]]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in lines[2:2 + n]])
    rcov = {"H": .31, "B": .84, "C": .76, "N": .71, "O": .66, "F": .57, "P": 1.07,
            "S": 1.05, "Cl": 1.02, "Ti": 1.6, "Br": 1.2, "I": 1.39}
    r = np.array([rcov.get(s, .77) for s in sym])
    metals = {"Ti", "Fe", "Zn", "Mg", "Na", "K"}
    adj = [[] for _ in range(n)]
    for i in range(n):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        cap = np.array([2.7 if (sym[i] in metals or s in metals) else 1.9 for s in sym])
        for j in np.where((d > 0.4) & (d < np.minimum(1.15 * (r + r[i]), cap)))[0]:
            if j != i:
                adj[i].append(int(j))
    seen = [False] * n
    comps = []
    for i in range(n):
        if seen[i]:
            continue
        stack, cur = [i], []
        while stack:
            k = stack.pop()
            if seen[k]:
                continue
            seen[k] = True
            cur.append(k)
            stack += adj[k]
        comps.append(cur)
    comps.sort(key=len, reverse=True)
    carrier, drug = comps[0], (comps[1] if len(comps) > 1 else [])
    d = out_dir or tempfile.mkdtemp()
    paths = []
    for name, idx in (("carrier", carrier), ("drug", drug)):
        p = os.path.join(d, f"_{name}.xyz")
        with open(p, "w") as fh:
            fh.write(f"{len(idx)}\n{name}\n")
            for k in idx:
                fh.write(f"{sym[k]:2s} {xyz[k,0]:14.8f} {xyz[k,1]:14.8f} {xyz[k,2]:14.8f}\n")
        paths.append(p)
    return paths[0], paths[1]


# --------------------------------------------------------------- static figures
def pocket_figure(receptor, ligand, out_png, key_res=None, size=(1600, 1300),
                  lig_carbon="limon", cartoon_color="skyblue", surface=True,
                  surf_transp=0.72, zoom_buffer=8.0, dpi=300, ray=1, turn=(-14, 5)):
    """Journal-style pocket figure: soft cartoon (+ optional see-through pocket
    surface), ligand as thick coloured sticks, key residues as thin grey sticks,
    polar contacts as yellow dashes. Residue names go in the figure caption,
    not on the render."""
    kr = key_res or []
    sel_res = " or ".join(f"resi {r}" for r in kr) if kr else "none"
    body = f"""
    cmd.load({_q(receptor)}, 'rec')
    cmd.load({_q(ligand)}, 'lig')
    cmd.remove('solvent')
    cmd.hide('everything')
    cmd.show('cartoon', 'rec and polymer')
    cmd.color({cartoon_color!r}, 'rec and polymer')
    cmd.show('sticks', 'lig')
    util.cnc('lig')
    cmd.color({lig_carbon!r}, 'lig and elem C')
    cmd.set('stick_radius', 0.26, 'lig')
    if {sel_res!r} != 'none':
        cmd.show('sticks', f'byres (rec and ({sel_res})) and (sidechain or name CA)')
        util.cnc(f'rec and ({sel_res})')
        cmd.color('grey60', f'rec and ({sel_res}) and elem C')
        cmd.set('stick_radius', 0.12, f'rec and ({sel_res})')
        cmd.distance('pol', 'lig', f'(rec and ({sel_res}) and (elem N+O))', 3.6, mode=2)
        cmd.hide('labels', 'pol')
    if {bool(surface)}:
        cmd.set('surface_mode', 0)
        cmd.set('surface_cavity_mode', 0)
        cmd.show('surface', 'rec and polymer')
        cmd.set('transparency', {surf_transp})
        cmd.set('surface_color', 'grey80')
    cmd.orient('lig')
    cmd.zoom('(rec and polymer within 13 of lig) or lig', {zoom_buffer})
    cmd.turn('y', {turn[0]}); cmd.turn('x', {turn[1]})
    cmd.png({_q(out_png)}, width={size[0]}, height={size[1]}, dpi={dpi}, ray={ray})
    """
    _run(body)
    return out_png


def complex_figure(struct, out_png, size=(1500, 1200), carbon="grey50",
                   mode="ball_stick", dpi=300, ray=1, tilt=18, turn=(0, 0, 0),
                   buffer=2.6, orient_on_carrier=True):
    """A drug-carrier / cluster complex as ray-traced ball-and-stick.

    orient_on_carrier: frame the view on the LARGER connected fragment (the 2D
    carrier), so the sheet sits roughly in-plane with the drug adsorbed on top,
    rather than edge-on. `tilt` then rocks it into a 3/4 view."""
    if mode == "cpk":
        show = "cmd.show('spheres')\ncmd.set('sphere_scale', 0.62)"
    elif mode == "sticks":
        show = "cmd.show('sticks')"
    else:
        show = "cmd.show('sticks')\ncmd.show('spheres')"
    orient_block = (
        "cmd.select('_f1', 'bymolecule (m and index 1)')\n"
        "cmd.select('_f2', 'm and not _f1')\n"
        "big = '_f1' if cmd.count_atoms('_f1') >= cmd.count_atoms('_f2') else '_f2'\n"
        "cmd.orient(big + ' and not hydro') if cmd.count_atoms('_f2') else cmd.orient('m and not hydro')\n"
        if orient_on_carrier else "cmd.orient('m and not hydro')\n")
    body = (
        f"cmd.load({_q(struct)}, 'm')\ncmd.hide('everything')\n"
        # PyMOL under-bonds metal-ligand contacts in bare .xyz files -> add them
        "for _mt, _cut in (('Ti', 2.35), ('Fe', 2.4), ('Zn', 2.4), ('Mg', 2.5)):\n"
        "    if cmd.count_atoms(f'm and elem {_mt}'):\n"
        "        cmd.bond(f'm and elem {_mt}', f'm and (elem O+C+N) within {_cut} of (m and elem {_mt})')\n"
        "util.cbaw('m')\n"
        f"cmd.color({carbon!r}, 'm and elem C')\n"
        "cmd.set('sphere_scale', 0.22)\ncmd.set('stick_radius', 0.135)\n"
        "cmd.set('valence', 1)\ncmd.set('valence_mode', 0)\n"
        + show + "\n"
        + orient_block +
        f"cmd.turn('x', {tilt})\n"
        f"cmd.turn('x', {turn[0]})\ncmd.turn('y', {turn[1]})\ncmd.turn('z', {turn[2]})\n"
        f"cmd.zoom('m', {buffer})\n"
        f"cmd.png({_q(out_png)}, width={size[0]}, height={size[1]}, dpi={dpi}, ray={ray})\n")
    _run(body)
    return out_png


def superpose_figure(receptor, pose_a, pose_b, out_png, size=(1600, 1300),
                     ca="palegreen", cb="lightorange", label_a="crystallographic",
                     label_b="redocked", surface=True, dpi=300):
    body = f"""
    cmd.load({_q(receptor)}, 'rec')
    cmd.load({_q(pose_a)}, 'pa')
    cmd.load({_q(pose_b)}, 'pb')
    cmd.remove('solvent')
    cmd.hide('everything')
    cmd.show('cartoon', 'rec and polymer'); cmd.color('grey80', 'rec')
    cmd.set('cartoon_transparency', 0.35, 'rec')
    for obj, c in (('pa', {ca!r}), ('pb', {cb!r})):
        cmd.show('sticks', obj); util.cnc(obj); cmd.color(c, obj + ' and elem C')
        cmd.set('stick_radius', 0.18, obj)
    if {bool(surface)}:
        cmd.show('surface', 'rec and polymer')
        cmd.set('transparency', 0.6, 'rec'); cmd.set('surface_color', 'grey90')
    cmd.orient('pa or pb'); cmd.zoom('pa or pb', 4.0)
    cmd.turn('y', -15)
    cmd.png({_q(out_png)}, width={size[0]}, height={size[1]}, dpi={dpi}, ray=1)
    """
    _run(body)
    return out_png


# ------------------------------------------------------------------- movies
def _encode(frames_dir, out_mp4, fps=30):
    pat = os.path.join(frames_dir, "frame%04d.png")
    subprocess.run([FFMPEG, "-y", "-framerate", str(fps), "-i", pat,
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p",
                    "-c:v", "libx264", "-crf", "18", "-preset", "slow", out_mp4],
                   capture_output=True, text=True, check=True)
    return out_mp4


def turntable(struct, out_mp4, n=150, fps=30, size=(880, 880), style="auto",
              carbon="grey55", ligand=None, ray=0):
    d = tempfile.mkdtemp()
    is_prot = str(struct).lower().endswith((".pdb", ".cif")) and style != "smallmol"
    if is_prot:
        setup = (f"cmd.load({_q(struct)}, 'm')\ncmd.remove('solvent')\ncmd.hide('everything')\n"
                 "cmd.show('cartoon', 'polymer')\n"
                 "cmd.spectrum('count', 'rainbow', 'polymer & name CA')\n"
                 "cmd.show('sticks', 'not polymer')\nutil.cnc('not polymer')\n"
                 "cmd.set('stick_radius', 0.22, 'not polymer')")
    else:
        setup = (f"cmd.load({_q(struct)}, 'm')\ncmd.hide('everything')\n"
                 "cmd.show('sticks')\ncmd.show('spheres')\nutil.cbaw('m')\n"
                 f"cmd.color({carbon!r}, 'm and elem C')\n"
                 "cmd.set('sphere_scale', 0.22)\ncmd.set('stick_radius', 0.135)")
    if ligand:
        setup += (f"\ncmd.load({_q(ligand)}, 'lg')\ncmd.show('sticks','lg')\n"
                  "cmd.show('spheres','lg')\nutil.cbay('lg')")
    body = (setup + "\n"
            "cmd.orient('m and not hydro')\ncmd.turn('x', 12)\ncmd.zoom('all', 2.2)\n"
            f"for i in range({n}):\n"
            f"    cmd.turn('y', 360.0 / {n})\n"
            f"    cmd.png(os.path.join({_q(d)}, 'frame%04d.png' % i), "
            f"width={size[0]}, height={size[1]}, dpi=150, ray={ray})\n")
    _run(body, timeout=3000)
    _encode(d, out_mp4, fps)
    shutil.rmtree(d, ignore_errors=True)
    return out_mp4


def _lift_vector(carrier_xyz_path, drug_xyz_path):
    """Direction to pull the drug away from the carrier: the carrier's plane
    normal (or centroid-difference for a non-planar carrier), signed toward the
    drug."""
    import numpy as np

    def rd(p):
        L = open(p).read().splitlines()
        n = int(L[0].split()[0])
        return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    cx, dx = rd(carrier_xyz_path), rd(drug_xyz_path)
    cc = cx.mean(0)
    _, s, vt = np.linalg.svd(cx - cc, full_matrices=False)
    normal = vt[2]
    planar = s[2] / (s[0] + 1e-9) < 0.35
    v = normal if planar else (dx.mean(0) - cc)
    if np.dot(v, dx.mean(0) - cc) < 0:
        v = -v
    n = np.linalg.norm(v)
    return (v / n) if n > 1e-6 else np.array([0.0, 0.0, 1.0])


def approach_movie(carrier, mobile, out_mp4, n_in=64, n_hold=14, n_spin=110,
                   fps=30, size=(880, 880), approach=15.0, carbon="grey55",
                   mob_carbon="limon", carrier_is_protein=False, ray=0):
    """`mobile` starts `approach` A off its bound pose along the carrier plane
    normal and eases in; then the assembled complex spins. The bound geometry is
    real - only the pre-binding offset is a synthetic 'approach'."""
    d = tempfile.mkdtemp()
    lv = _lift_vector(carrier, mobile)
    if carrier_is_protein:
        carr = ("cmd.show('cartoon','polymer')\ncmd.color('skyblue','polymer')\n"
                "cmd.set('cartoon_transparency',0.12)\n"
                "cmd.show('surface','polymer')\ncmd.set('transparency',0.78)\n"
                "cmd.set('surface_color','grey80')")
    else:
        carr = ("cmd.show('sticks','car')\ncmd.show('spheres','car')\nutil.cbaw('car')\n"
                f"cmd.color({carbon!r},'car and elem C')")
    body = (
        f"cmd.load({_q(carrier)}, 'car')\ncmd.load({_q(mobile)}, 'mob')\n"
        "cmd.remove('solvent')\ncmd.hide('everything')\n" + carr + "\n"
        "cmd.show('sticks','mob')\ncmd.show('spheres','mob')\nutil.cbaw('mob')\n"
        f"cmd.color({mob_carbon!r}, 'mob and elem C')\n"
        "cmd.set('sphere_scale', 0.22)\ncmd.set('stick_radius', 0.145)\n"
        "import numpy as _np\n"
        f"v = _np.array([{float(lv[0])!r}, {float(lv[1])!r}, {float(lv[2])!r}])\n"
        "mc = _np.array(cmd.centerofmass('mob'))\n"
        # frame on the carrier so the sheet sits in-plane; keep this view fixed
        "cmd.orient('car and not hydro')\ncmd.turn('x', 14)\n"
        # a pseudoatom at the fully-displaced drug position so zoom fits the whole path
        f"cmd.pseudoatom('_far', pos=list(mc + v*{approach}))\n"
        f"cmd.zoom('car or mob or _far', 1.5)\ncmd.delete('_far')\n"
        f"NI, NS, APP = {n_in}, {n_spin}, {approach}\n"
        "def ease(t): return t*t*(3-2*t)\n"
        "cmd.translate(list(v*APP), 'mob', camera=0)\n"
        "prev = 1.0\nfi = 0\n"
        "for i in range(NI):\n"
        "    cur = 1.0 - ease(i/(NI-1))\n"
        "    cmd.translate(list(v*APP*(cur-prev)), 'mob', camera=0)\n"
        "    prev = cur\n"
        f"    cmd.png(os.path.join({_q(d)}, 'frame%04d.png' % fi), width={size[0]}, height={size[1]}, dpi=150, ray={ray})\n"
        "    fi += 1\n"
        f"for _ in range({n_hold}):\n"
        f"    cmd.png(os.path.join({_q(d)}, 'frame%04d.png' % fi), width={size[0]}, height={size[1]}, dpi=150, ray={ray})\n"
        "    fi += 1\n"
        "for i in range(NS):\n"
        "    cmd.turn('y', 360.0 / NS)\n"
        f"    cmd.png(os.path.join({_q(d)}, 'frame%04d.png' % fi), width={size[0]}, height={size[1]}, dpi=150, ray={ray})\n"
        "    fi += 1\n")
    _run(body, timeout=3000)
    _encode(d, out_mp4, fps)
    shutil.rmtree(d, ignore_errors=True)
    return out_mp4

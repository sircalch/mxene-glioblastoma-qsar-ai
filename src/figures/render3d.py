"""
render3d.py - ray-traced molecular renders (open-source PyMOL, headless).

PyMOL lives in a separate conda env (PYMOL_PY); each render is a small script
run in a subprocess. Residue labels are NOT drawn by PyMOL: the script
projects the label anchor atoms to image pixels and returns them as JSON, so
labels are typeset by matplotlib in the figure font.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import textwrap

PYMOL_PY = os.environ.get("PYMOL_PY", r"C:\Users\Andre\mm\mol\python.exe")
AVAILABLE = os.path.exists(PYMOL_PY)

_HEADER = r"""
import os, sys, json
import numpy as np
os.environ.setdefault('PYMOL_PATH', r'{root}\share\pymol')
import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd, util
cmd.bg_color('white')
cmd.set('ray_opaque_background', 0)
cmd.set('orthoscopic', 1)
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_gain', 0.05)
cmd.set('ray_trace_color', 'grey30')
cmd.set('ray_shadows', 1)
cmd.set('ambient', 0.35)
cmd.set('direct', 0.55)
cmd.set('reflect', 0.25)
cmd.set('specular', 0.25)
cmd.set('spec_power', 120)
cmd.set('light_count', 2)
cmd.set('depth_cue', 1)
cmd.set('fog_start', 0.55)
cmd.set('antialias', 2)
cmd.set('cartoon_fancy_helices', 1)
cmd.set('cartoon_smooth_loops', 1)
cmd.set('cartoon_highlight_color', -1)
cmd.set('cartoon_side_chain_helper', 1)
cmd.set('surface_quality', 1)
cmd.set('two_sided_lighting', 1)
cmd.set('valence', 0)
cmd.set('dash_gap', 0.35)
cmd.set('dash_radius', 0.07)
cmd.set('dash_color', 'grey20')

def face_camera(target, context):
    '''Rotate the view so the vector context-centre -> target-centre points at the viewer.'''
    import math
    for _ in range(3):
        v = cmd.get_view()
        R = np.array(v[0:9]).reshape(3, 3)
        d = np.array(cmd.get_coords(target)).mean(0) - np.array(cmd.get_coords(context)).mean(0)
        c = d @ R
        cmd.turn('y', -math.degrees(math.atan2(c[0], c[2])))
        v = cmd.get_view()
        R = np.array(v[0:9]).reshape(3, 3)
        c = d @ R
        cmd.turn('x', math.degrees(math.atan2(c[1], c[2])))


def project(sel, W, H):
    '''Pixel coordinates of the centroid of `sel` in a W x H orthoscopic render.'''
    v = cmd.get_view()
    R = np.array(v[0:9]).reshape(3, 3)
    cam = np.array(v[9:12]); org = np.array(v[12:15])
    xyz = np.array(cmd.get_coords(sel)).mean(0)
    c = (xyz - org) @ R + cam
    half_h = -cam[2] * np.tan(np.radians(float(cmd.get('field_of_view'))) / 2)
    half_w = half_h * W / H
    return [float(W / 2 + c[0] / half_w * W / 2), float(H / 2 - c[1] / half_h * H / 2)]
""".replace("{root}", os.path.dirname(os.path.dirname(PYMOL_PY)))


def _run(body: str, timeout: int = 1800) -> dict:
    """Run a PyMOL script; whatever it passes to emit() comes back as a dict
    (PyMOL swallows stdout, so results travel through a temp JSON file)."""
    fd, path = tempfile.mkstemp(suffix=".py")
    res = path[:-3] + ".json"
    head = _HEADER + f"\ndef emit(d):\n    json.dump(d, open({res.replace(os.sep, '/')!r}, 'w'))\n"
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(head + textwrap.dedent(body))
    try:
        p = subprocess.run([PYMOL_PY, path], capture_output=True, text=True, timeout=timeout)
        if p.returncode != 0 or "Traceback" in p.stderr:
            raise RuntimeError(f"PyMOL failed:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
        if os.path.exists(res):
            with open(res) as fh:
                return json.load(fh)
        return {}
    finally:
        os.unlink(path)
        if os.path.exists(res):
            os.unlink(res)


def _q(p: str) -> str:
    return repr(os.path.abspath(p).replace("\\", "/"))


def first_model_pdb(pdbqt: str, out_pdb: str) -> str:
    """MODEL 1 of a Vina .pdbqt -> plain PDB (element column restored)."""
    lines = []
    with open(pdbqt, encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith("ENDMDL"):
                break
            if ln.startswith(("ATOM", "HETATM")):
                el = ln[77:79].strip() or ln[12:14].strip()
                el = {"A": "C", "OA": "O", "NA": "N", "SA": "S", "HD": "H"}.get(el, el)
                lines.append("HETATM" + ln[6:17] + "LIG L 999    " + ln[30:54]
                             + "  1.00 20.00          " + f"{el:>2}\n")
    with open(out_pdb, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
        fh.write("END\n")
    return out_pdb


def protein_overview(receptor: str, ligand: str, out_png: str, size=(1800, 1500)):
    """Whole receptor: soft cartoon per chain + faint surface, ligand as CPK."""
    W, H = size
    body = f"""
    cmd.load({_q(receptor)}, 'rec'); cmd.load({_q(ligand)}, 'lig')
    cmd.remove('solvent or inorganic or (rec and not polymer)')
    cmd.hide('everything')
    cmd.show('cartoon', 'rec')
    cmd.set_color('chA', [0.62, 0.72, 0.84]); cmd.set_color('chB', [0.84, 0.86, 0.89])
    cmd.color('chA', 'rec and chain A'); cmd.color('chB', 'rec and chain B')
    cmd.show('surface', 'rec')
    cmd.set('transparency', 0.82)
    cmd.set('surface_color', 'grey90')
    cmd.show('spheres', 'lig and not elem H')
    util.cnc('lig'); cmd.set_color('ligC', [0.91, 0.55, 0.16]); cmd.color('ligC', 'lig and elem C')
    cmd.orient('rec')
    cmd.turn('x', -12)
    cmd.zoom('rec', 1.0)
    cmd.png({_q(out_png)}, width={W}, height={H}, dpi=600, ray=1)
    emit(({{'lig': project('lig', {W}, {H})}}))
    """
    return _run(body)


def pocket_closeup(receptor, ligands, out_png, cutoff=4.0, size=(1800, 1500), turn=(0, 0),
                   cofactors=None, chain="A", label_min_contacts=1, zoom=1.5, slab=40):
    """Binding-site close-up.

    ligands: list of dicts, first one defines the site and the polar contacts:
        {"path": file | "sel": selection inside `receptor`, "color": (r,g,b),
         "radius": stick radius, "name": object name}
    cofactors: selection inside `receptor` shown as thin sticks (+ ions as spheres).
    Returns {'residues': {label: [px, py]}, 'ligands': {name: [px, py]}}."""
    W, H = size
    lig_lines = []
    for i, L in enumerate(ligands):
        nm = L.get("name", f"lig{i}")
        if "path" in L:
            lig_lines.append(f"cmd.load({_q(L['path'])}, '{nm}')")
        else:
            lig_lines.append(f"cmd.create('{nm}', 'rec and chain {chain} and ({L['sel']})')")
            lig_lines.append(f"cmd.remove('rec and ({L['sel']})')")
        c = L.get("color", (0.91, 0.55, 0.16))
        lig_lines += [f"cmd.show('sticks', '{nm} and not hydro')", f"util.cnc('{nm}')",
                      f"cmd.set_color('c_{nm}', {list(c)})", f"cmd.color('c_{nm}', '{nm} and elem C')",
                      f"cmd.set('stick_radius', {L.get('radius', 0.24)}, '{nm}')"]
    names = [L.get("name", f"lig{i}") for i, L in enumerate(ligands)]
    first = names[0]
    cof = cofactors or "none"
    body = f"""
    cmd.load({_q(receptor)}, 'rec')
    cmd.remove('solvent or (rec and not chain {chain})')
    cmd.remove('hydro and not (elem H and neighbor elem N+O)')
    cmd.hide('everything')
    """ + "\n    ".join(lig_lines) + f"""
    cmd.select('cof', 'rec and ({cof})')
    cmd.select('site', 'byres (rec and polymer and not cof) within {cutoff} of {first}')
    cmd.show('cartoon', 'rec and polymer')
    cmd.set('cartoon_transparency', 0.72)
    cmd.set_color('cart', [0.78, 0.83, 0.89]); cmd.color('cart', 'rec and polymer')
    cmd.show('sticks', 'site and not name N+C+O')
    util.cnc('site'); cmd.set_color('resC', [0.62, 0.66, 0.72]); cmd.color('resC', 'site and elem C')
    cmd.set('stick_radius', 0.14, 'site')
    if cmd.count_atoms('cof'):
        cmd.show('sticks', 'cof and not elem Mg+Zn+Na+Ca')
        util.cnc('cof'); cmd.set_color('cofC', [0.55, 0.70, 0.55]); cmd.color('cofC', 'cof and elem C')
        cmd.set('stick_radius', 0.12, 'cof')
        cmd.show('spheres', 'cof and elem Mg+Zn+Na+Ca'); cmd.set('sphere_scale', 0.35, 'cof')
        cmd.set_color('ion', [0.45, 0.80, 0.45]); cmd.color('ion', 'cof and elem Mg+Zn+Na+Ca')
    cmd.distance('hb', '{first} and (elem N+O)', '(site or cof) and (elem N+O)', 3.5, mode=2)
    cmd.hide('labels', 'hb')
    cmd.show('surface', 'site')
    cmd.set('transparency', 0.80)
    cmd.set('surface_color', 'grey85')
    cmd.orient('{first}')
    face_camera('{first}', 'rec and polymer')
    cmd.turn('y', {turn[0]}); cmd.turn('x', {turn[1]})
    cmd.zoom('{first} or site', {zoom})
    cmd.clip('slab', {slab})
    cmd.png({_q(out_png)}, width={W}, height={H}, dpi=600, ray=1)
    labs = {{}}
    stored = []
    cmd.iterate('site and name CA', 'stored.append((chain, resn, resi))', space={{'stored': stored}})
    for ch, rn, ri in stored:
        sel = f'site and chain {{ch}} and resi {{ri}} and not name N+C+O+CA'
        if cmd.count_atoms(sel) == 0:
            sel = f'site and chain {{ch}} and resi {{ri}} and name CA'
        labs[f'{{rn.capitalize()}}{{ri}}'] = project(sel, {W}, {H})
    emit({{'residues': labs, 'ligands': {{n: project(n, {W}, {H}) for n in {names!r}}}}})
    """
    return _run(body)


def molecule(struct: str, out_png: str, size=(1400, 1150), carbon=(0.42, 0.45, 0.50),
             tilt=18, turn=(0, 0, 0), buffer=2.6, orient_on_carrier=True):
    """Ball-and-stick render of a drug/cage complex or an isolated cage."""
    body = f"""
    cmd.load({_q(struct)}, 'm')
    cmd.hide('everything')
    util.cbaw('m')
    cmd.set_color('cC', {list(carbon)}); cmd.color('cC', 'm and elem C')
    cmd.set_color('cB', [0.96, 0.71, 0.66]); cmd.color('cB', 'm and elem B')
    cmd.set_color('cN', [0.20, 0.33, 0.85]); cmd.color('cN', 'm and elem N')
    cmd.set_color('cTi', [0.62, 0.68, 0.80]); cmd.color('cTi', 'm and elem Ti')
    cmd.set('sphere_scale', 0.22); cmd.set('stick_radius', 0.13)
    cmd.show('sticks'); cmd.show('spheres')
    cmd.set('depth_cue', 0)
    cmd.select('_f1', 'bymolecule (m and index 1)')
    cmd.select('_f2', 'm and not _f1')
    big = '_f1' if cmd.count_atoms('_f1') >= cmd.count_atoms('_f2') else '_f2'
    cmd.orient((big + ' and not hydro') if ({orient_on_carrier} and cmd.count_atoms('_f2')) else 'm and not hydro')
    cmd.turn('x', {tilt})
    cmd.turn('x', {turn[0]}); cmd.turn('y', {turn[1]}); cmd.turn('z', {turn[2]})
    cmd.zoom('m', {buffer})
    cmd.png({_q(out_png)}, width={size[0]}, height={size[1]}, dpi=600, ray=1)
    """
    _run(body)
    return out_png

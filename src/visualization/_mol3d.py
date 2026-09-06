"""
_mol3d.py - publication-quality 3D molecular rendering via PyVista (off-screen).

Reads .xyz / .pdb / .pdbqt, perceives bonds from covalent radii, and renders
ball-and-stick (small molecules / clusters) and translucent molecular surfaces
(proteins) with a 3-point light rig, physically-based materials and SSAA.

No GPU required: VTK's off-screen path works in this environment.
"""
from __future__ import annotations
import os
import numpy as np

try:
    import pyvista as pv
    pv.OFF_SCREEN = True
except Exception as exc:  # pragma: no cover
    pv = None
    _IMPORT_ERROR = exc

# ---------------------------------------------------------------- element data
# Covalent radii (Angstrom, Cordero 2008, trimmed) and CPK-ish colours tuned
# for a bright white background.
COV = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Ti": 1.60, "Mg": 1.41, "Na": 1.66, "K": 2.03, "Fe": 1.52,
}
COL = {
    "H": "#EAEAEA", "B": "#EF9BB6", "C": "#4D4D4D", "N": "#2F62C4", "O": "#E23B30",
    "F": "#8FD14F", "P": "#F0912B", "S": "#E3B93B", "Cl": "#49B849", "Br": "#A0522D",
    "I": "#7A2FA0", "Ti": "#9AA0A6", "Mg": "#3FBF6F", "Na": "#8A6FE8", "K": "#7A52D0",
    "Fe": "#D77A2B", "X": "#B0B0B0",
}
# AutoDock PDBQT atom types -> element
_ADT = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O",
        "SA": "S", "S": "S", "HD": "H", "HS": "H", "H": "H", "F": "F", "Cl": "Cl",
        "CL": "Cl", "Br": "Br", "BR": "Br", "I": "I", "P": "P", "Mg": "Mg",
        "MG": "Mg", "Zn": "Zn", "ZN": "Zn", "Ca": "Ca", "CA": "Ca", "Fe": "Fe",
        "FE": "Fe", "MN": "Mn", "e": "C", "G0": "C", "CG0": "C"}
VDW = {"H": 1.10, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
       "P": 1.80, "S": 1.80, "Cl": 1.75, "Ti": 2.15, "Mg": 1.73}


# --------------------------------------------------------------------- readers
def _read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    sym, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        sym.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return sym, np.asarray(xyz, float)


def _elem_from_pdb_line(ln):
    e = ln[76:78].strip()
    if e:
        return e[0].upper() + e[1:].lower()
    nm = ln[12:16].strip()
    nm = "".join(c for c in nm if c.isalpha())
    if len(nm) >= 2 and nm[:2].capitalize() in COV:
        return nm[:2].capitalize()
    return nm[:1].upper() if nm else "X"


def _read_pdb(path, hetatm=True, atom=True, chain=None, resnames=None,
              max_models=1):
    sym, xyz = [], []
    models = 0
    for ln in open(path):
        rec = ln[:6].strip()
        if rec == "MODEL":
            models += 1
            if models > max_models:
                break
        if rec == "ENDMDL":
            if models >= max_models:
                break
        if rec not in ("ATOM", "HETATM"):
            continue
        if rec == "ATOM" and not atom:
            continue
        if rec == "HETATM" and not hetatm:
            continue
        if rec == "HETATM" and ln[17:20].strip() in ("HOH", "WAT", "DOD"):
            continue
        if chain and ln[21] != chain:
            continue
        if resnames and ln[17:20].strip() not in resnames:
            continue
        try:
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        except ValueError:
            continue
        sym.append(_elem_from_pdb_line(ln))
        xyz.append([x, y, z])
    return sym, np.asarray(xyz, float) if xyz else np.zeros((0, 3))


def _read_pdbqt(path, model=1):
    sym, xyz = [], []
    cur = 0
    for ln in open(path):
        if ln.startswith("MODEL"):
            cur = int(ln.split()[1])
        if cur and cur != model:
            continue
        if ln[:6].strip() not in ("ATOM", "HETATM"):
            continue
        try:
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        except ValueError:
            continue
        adt = ln[76:79].strip() or ln[12:16].strip()
        e = _ADT.get(adt) or _ADT.get(adt.upper()) or _ADT.get(adt.capitalize())
        if not e:
            e = _elem_from_pdb_line(ln)
        sym.append(e if e in COV else "C")
        xyz.append([x, y, z])
    return sym, np.asarray(xyz, float)


def load(path, **kw):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".xyz":
        return _read_xyz(path)
    if ext == ".pdbqt":
        return _read_pdbqt(path, model=kw.get("model", 1))
    return _read_pdb(path, **{k: v for k, v in kw.items()
                              if k in ("hetatm", "atom", "chain",
                                       "resnames", "max_models")})


# ----------------------------------------------------------------- geometry
_METALS = {"Ti", "Fe", "Zn", "Mg", "Na", "K", "Ca", "Mn", "Cu", "Co", "Ni"}


def bonds(sym, xyz, scale=1.15):
    n = len(sym)
    r = np.array([COV.get(s, 0.77) for s in sym])
    out = []
    for i in range(n):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        metal_i = sym[i] in _METALS
        cap = np.array([2.75 if (metal_i or s in _METALS) else 1.95 for s in sym])
        cut = np.minimum(scale * (r + r[i]), cap)
        for j in np.where((d > 0.4) & (d < cut))[0]:
            j = int(j)
            if j <= i:
                continue
            if "H" in (sym[i], sym[j]) and d[j] > 1.30:
                continue
            out.append((i, j))
    return out


# ----------------------------------------------------------------- rendering
def _light_rig(pl):
    pl.remove_all_lights()
    # soft key + two fills + gentle rim; all headlight-relative so the molecule
    # is lit from the camera's frame no matter the auto-orientation.
    for pos, inten, col in [
        ((-0.6, 0.8, 1.0), 1.05, "white"),
        ((0.9, 0.2, 0.5), 0.42, "#e9f0ff"),
        ((-0.3, -0.9, 0.2), 0.30, "#fff4e8"),
        ((0.0, 0.1, -1.0), 0.22, "white"),
    ]:
        lt = pv.Light(position=pos, light_type="camera light",
                      intensity=inten, color=col)
        pl.add_light(lt)


def _atom_radius(s, ball):
    if s == "H":
        return ball * 0.55
    if s in ("Ti", "I", "Br", "K", "Na"):
        return ball * 1.5
    if s in ("P", "S", "Cl", "B"):
        return ball * 1.15
    return ball


def _add_ball_stick(pl, sym, xyz, ball=0.30, stick=0.105, bond_scale=1.18,
                    carbon=None, opacity=1.0):
    if len(sym) == 0:
        return
    carbon = carbon or COL["C"]
    b = bonds(sym, xyz, bond_scale)
    # bonds: one tube per half-bond, coloured by the nearer atom
    for i, j in b:
        mid = (xyz[i] + xyz[j]) / 2
        for a, m in ((i, mid), (j, mid)):
            v = m - xyz[a]
            h = float(np.linalg.norm(v))
            if h < 1e-3:
                continue
            cyl = pv.Cylinder(center=tuple((xyz[a] + m) / 2), direction=tuple(v),
                              radius=stick, height=h, resolution=20, capping=True)
            c = carbon if sym[a] == "C" else COL.get(sym[a], COL["X"])
            pl.add_mesh(cyl, color=c, pbr=True, metallic=0.0, roughness=0.5,
                        smooth_shading=True, opacity=opacity, specular=0.15)
    # atoms, glyphed per element for speed
    unit = pv.Sphere(radius=1.0, theta_resolution=44, phi_resolution=44)
    for s in sorted(set(sym)):
        idx = [k for k in range(len(sym)) if sym[k] == s]
        pd = pv.PolyData(xyz[idx])
        glyph = pd.glyph(geom=unit, scale=False, factor=_atom_radius(s, ball),
                         orient=False)
        c = carbon if s == "C" else COL.get(s, COL["X"])
        pl.add_mesh(glyph, color=c, pbr=True, metallic=0.0,
                    roughness=0.34 if s != "H" else 0.55,
                    smooth_shading=True, opacity=opacity, specular=0.25,
                    specular_power=18)


def _add_surface(pl, sym, xyz, color="#c9d3e0", opacity=0.42,
                 probe=1.3, grid=0.45, iso=2.4):
    """Smooth Gaussian-blur molecular surface (Blinn-style), heavily relaxed so
    it reads as a clean SES rather than a lumpy blob."""
    if len(sym) < 4:
        return
    r = np.array([VDW.get(s, 1.7) for s in sym])
    pad = r.max() + probe + 2.5
    lo = xyz.min(0) - pad
    hi = xyz.max(0) + pad
    dims = np.maximum(((hi - lo) / grid).astype(int), 12)
    nx, ny, nz = dims
    gx = np.linspace(lo[0], hi[0], nx)
    gy = np.linspace(lo[1], hi[1], ny)
    gz = np.linspace(lo[2], hi[2], nz)
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    dens = np.zeros_like(X)
    for k in range(len(sym)):
        s = r[k] + probe
        d2 = (X - xyz[k, 0])**2 + (Y - xyz[k, 1])**2 + (Z - xyz[k, 2])**2
        dens += np.exp(-1.8 * d2 / (s * s))
    g = pv.ImageData(dimensions=(nx, ny, nz),
                     spacing=(gx[1]-gx[0], gy[1]-gy[0], gz[1]-gz[0]),
                     origin=tuple(lo))
    g["d"] = dens.flatten(order="F")
    try:
        surf = g.contour([iso], scalars="d")
        surf = surf.smooth_taubin(n_iter=60, pass_band=0.05)
        surf = surf.compute_normals(auto_orient_normals=True, split_vertices=False)
        pl.add_mesh(surf, color=color, opacity=opacity, pbr=True, metallic=0.0,
                    roughness=0.6, smooth_shading=True, specular=0.12,
                    diffuse=0.9)
    except Exception:
        pass


# ------------------------------------------------------ protein cartoon ribbon
def _pdb_with_ss(pdb_path):
    """Rewrite a PDB with proper HELIX/SHEET records (secondary structure from
    biotite's P-SEA) so vtkPDBReader/vtkProteinRibbonFilter draw real
    helices and arrows. Returns (tmp_path, ca_xyz)."""
    import tempfile
    import biotite.structure.io.pdb as _pdb
    import biotite.structure as _st
    arr = _pdb.PDBFile.read(pdb_path).get_structure(model=1)
    arr = arr[_st.filter_amino_acids(arr)]
    recs = []
    hn = sn = 0
    for ch in sorted(set(arr.chain_id)):
        c = arr[arr.chain_id == ch]
        try:
            sse = _st.annotate_sse(c)
        except Exception:
            sse = np.array(["c"] * len(np.unique(c.res_id)))
        res = np.unique(c.res_id)
        rn = {int(a.res_id): a.res_name for a in c[c.atom_name == "CA"]}
        i = 0
        while i < len(sse):
            j = i
            while j < len(sse) and sse[j] == sse[i]:
                j += 1
            a, b = int(res[i]), int(res[j-1])
            ra = rn.get(a, "ALA")[:3].rjust(3)
            rb = rn.get(b, "ALA")[:3].rjust(3)
            if sse[i] == "a" and j - i >= 4:
                hn += 1
                recs.append(f"HELIX  {hn:3d} {hn:3d} {ra} {ch} {a:4d}  {rb} {ch} {b:4d}  1"
                            .ljust(76))
            elif sse[i] == "b" and j - i >= 2:
                sn += 1
                recs.append(f"SHEET  {sn:3d} S{ch}{1:2d} {ra} {ch}{a:4d}  {rb} {ch}{b:4d}  0"
                            .ljust(70))
            i = j
    lines = _pdb.PDBFile()
    lines.set_structure(arr)
    atom_lines = [ln for ln in lines.lines if ln[:6] in ("ATOM  ", "HETATM", "TER   ")]
    tmp = os.path.join(tempfile.gettempdir(),
                       "_ss_" + os.path.basename(pdb_path))
    with open(tmp, "w") as fh:
        fh.write("\n".join(recs) + "\n")
        fh.write("\n".join(atom_lines) + "\nEND\n")
    ca = arr[arr.atom_name == "CA"].coord
    return tmp, np.asarray(ca, float)


def load_protein(pdb_path, chain=None):
    """Return a dict for add_cartoon: {'ss_pdb': path, 'ca': xyz}."""
    ss_pdb, ca = _pdb_with_ss(pdb_path)
    return {"ss_pdb": ss_pdb, "ca": ca, "src": pdb_path}


def _hex(c):
    c = c.lstrip("#")
    return np.array([int(c[i:i+2], 16) for i in (0, 2, 4)], dtype=np.uint8)


def add_cartoon(pl, prot, color_by="ss", helix="#E8776B", sheet="#E8C15A",
                loop="#F2F2F2"):
    """Secondary-structure cartoon via vtkProteinRibbonFilter (real helices,
    arrow strands, coil tube), recoloured to a clean journal palette.

    color_by: 'ss' (helix/sheet/loop palette) or a hex string (uniform)."""
    import vtk
    rd = vtk.vtkPDBReader()
    rd.SetFileName(prot["ss_pdb"])
    rd.Update()
    rib = vtk.vtkProteinRibbonFilter()
    rib.SetInputConnection(rd.GetOutputPort())
    try:
        rib.SetDrawSmallMoleculesAsSpheres(False)
    except Exception:
        pass
    rib.Update()
    mesh = pv.wrap(rib.GetOutput())
    if mesh.n_points == 0:
        return
    rgb = np.asarray(mesh.point_data.get("RGB"))
    if rgb is None or color_by != "ss":
        col = _hex(color_by) if isinstance(color_by, str) and color_by.startswith("#") else _hex(helix)
        pl.add_mesh(mesh, color=(col/255.0), pbr=True, metallic=0.0,
                    roughness=0.42, smooth_shading=True, specular=0.25,
                    specular_power=15)
        return
    # vtkProteinRibbonFilter marks helix ~ magenta/red, sheet ~ yellow, loop ~ white/grey.
    out = rgb.copy()
    r, g, bl = rgb[:, 0].astype(int), rgb[:, 1].astype(int), rgb[:, 2].astype(int)
    is_loop = (r > 190) & (g > 190) & (bl > 190)
    is_sheet = (r > 150) & (g > 120) & (bl < 130) & ~is_loop
    is_helix = ~is_loop & ~is_sheet
    out[is_helix] = _hex(helix)
    out[is_sheet] = _hex(sheet)
    out[is_loop] = _hex(loop)
    mesh.point_data["cartoon_rgb"] = out
    pl.add_mesh(mesh, scalars="cartoon_rgb", rgb=True, pbr=True, metallic=0.0,
                roughness=0.42, smooth_shading=True, specular=0.28,
                specular_power=15)



def _auto_view(all_xyz, view="face", roll=0.0):
    """Return [position, focal, viewup].

    view: "face" look down the thinnest axis (planar sheets / rings flat-on)
          "edge" look along a long axis, thinnest axis vertical (shows the
                 pi-pi stacking gap of a drug on a sheet)
          "3q"   three-quarter perspective
    """
    c = all_xyz.mean(0)
    X = all_xyz - c
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    e1, e2, e3 = vt[0], vt[1], vt[2]        # major, minor, thinnest
    span = float(np.linalg.norm(all_xyz.max(0) - all_xyz.min(0)))
    if view == "face":
        view_dir, up = e3, e2
    elif view == "edge":
        view_dir, up = e1, e3
    else:  # 3q
        view_dir = 0.6 * e3 + 0.5 * e1 - 0.35 * e2
        view_dir /= np.linalg.norm(view_dir)
        up = e3
    if roll:
        k = view_dir
        up = (up * np.cos(roll) + np.cross(k, up) * np.sin(roll)
              + k * np.dot(k, up) * (1 - np.cos(roll)))
    pos = c + view_dir * span * 2.4
    return [tuple(pos), tuple(c), tuple(up)]


def _autocrop(img, pad=12, bg=250):
    mask = np.any(img < bg, axis=2)
    if not mask.any():
        return img
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def _draw_spec(pl, sp, offset=None):
    """Render one spec dict into a plotter. Returns its xyz (for framing)."""
    off = np.zeros(3) if offset is None else np.asarray(offset, float)
    got = []
    if sp.get("cartoon"):
        add_cartoon(pl, sp["cartoon"], color_by=sp.get("cartoon_color", "ss"),
                    tube_r=sp.get("tube_r", 0.28))
        for seg in sp["cartoon"]["segments"]:
            got.append(seg["ca"])
    sym = sp.get("sym")
    if sym is not None and len(sym):
        xyz = np.asarray(sp["xyz"], float) + off
        got.append(xyz)
        if sp.get("surface"):
            _add_surface(pl, sym, xyz, color=sp.get("surf_color", "#c9d3e0"),
                         opacity=sp.get("surf_opacity", 0.42))
        if sp.get("style", "ball_stick") != "none":
            _add_ball_stick(pl, sym, xyz, ball=sp.get("ball", 0.30),
                            stick=sp.get("stick", 0.105),
                            bond_scale=sp.get("bond_scale", 1.12),
                            carbon=sp.get("carbon"),
                            opacity=sp.get("opacity", 1.0))
    return np.vstack(got) if got else np.zeros((0, 3))


def _finish(pl, allxyz, view, roll, zoom, ssao, ssaa, orthographic):
    _light_rig(pl)
    try:
        pl.enable_depth_peeling(12)
    except Exception:
        pass
    if ssao:
        try:
            pl.enable_ssao(radius=2.2, bias=0.5, kernel_size=128, blur=True)
        except Exception:
            pass
    if orthographic:
        pl.camera.SetParallelProjection(True)
    pl.camera_position = list(_auto_view(allxyz, view=view, roll=roll))
    pl.camera.zoom(zoom)
    if ssaa:
        try:
            pl.enable_anti_aliasing("ssaa")
        except Exception:
            pass


def render_array(specs, size=(1600, 1200), zoom=1.15, bg="white", ssaa=True,
                 orthographic=False, view="face", roll=0.0, ssao=True,
                 crop=True):
    """Render and return an (H, W, 3) uint8 RGB array.

    specs: list of dicts. Each may carry {"sym","xyz"} for ball&stick/surface
    (options: ball, stick, bond_scale, carbon, opacity, surface, surf_color,
    surf_opacity, style('none')) and/or {"cartoon": load_protein(...)} for a
    secondary-structure ribbon (cartoon_color: 'ss'|'chain'|'rainbow'|hex)."""
    if pv is None:  # pragma: no cover
        raise RuntimeError(f"pyvista unavailable: {_IMPORT_ERROR}")
    pl = pv.Plotter(off_screen=True, window_size=size, lighting="none")
    pl.set_background(bg)
    allxyz = [x for sp in specs for x in [_draw_spec(pl, sp)] if len(x)]
    allxyz = np.vstack(allxyz)
    _finish(pl, allxyz, view, roll, zoom, ssao, ssaa, orthographic)
    img = pl.screenshot(return_img=True)
    pl.close()
    if img.ndim == 3 and img.shape[2] == 4:
        img = img[:, :, :3]
    return _autocrop(img) if crop else img


def render(specs, out_path, **kw):
    import imageio.v2 as imageio
    img = render_array(specs, **kw)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    imageio.imwrite(out_path, img)
    return out_path


def _even(v):
    v = int(round(v))
    return v + (v & 1)


def _write_mp4(frames, out_mp4, fps):
    import imageio.v2 as imageio
    os.makedirs(os.path.dirname(os.path.abspath(out_mp4)) or ".", exist_ok=True)
    h, w = frames[0].shape[:2]
    tw, th = _even(w), _even(h)
    fr = [f[:th, :tw] if (f.shape[0] >= th and f.shape[1] >= tw) else f
          for f in frames]
    imageio.mimsave(out_mp4, fr, fps=fps, quality=9, codec="libx264",
                    macro_block_size=1,
                    output_params=["-pix_fmt", "yuv420p", "-crf", "18"])
    return out_mp4


def turntable(specs, out_mp4, size=(960, 960), n=140, zoom=1.25, fps=30,
              bg="white", view="3q", rock=False):
    """Rotating-structure movie (full 360 by default)."""
    if pv is None:  # pragma: no cover
        raise RuntimeError("pyvista unavailable")
    pl = pv.Plotter(off_screen=True, window_size=size, lighting="none")
    pl.set_background(bg)
    allxyz = np.vstack([_draw_spec(pl, sp) for sp in specs])
    _finish(pl, allxyz, view, 0.0, zoom, True, False, False)
    try:
        pl.enable_anti_aliasing("msaa")
    except Exception:
        pass
    step = 360.0 / n
    frames = []
    for i in range(n):
        ang = (step * (1 if not rock else np.sin(2*np.pi*i/n) * 0.6))
        pl.camera.Azimuth(ang)          # incremental VTK rotation
        pl.render()
        frames.append(np.asarray(pl.screenshot(return_img=True))[..., :3])
    pl.close()
    return _write_mp4(frames, out_mp4, fps)


def assembly_movie(scene_specs, mobile_idx, out_mp4, size=(960, 960),
                   approach=14.0, n_in=70, n_hold=18, n_spin=110, fps=30,
                   zoom=1.2, bg="white", view="3q", along=None):
    """Animate one spec (mobile_idx) sliding into its bound pose, then a spin.

    scene_specs[mobile_idx] must have {"sym","xyz"} at the BOUND geometry.
    The mobile fragment starts `approach` A away along `along` (default: the
    vector from the rest-of-scene centroid to the fragment centroid) and eases in.
    """
    if pv is None:  # pragma: no cover
        raise RuntimeError("pyvista unavailable")
    mob = scene_specs[mobile_idx]
    mob_xyz = np.asarray(mob["xyz"], float)
    rest = [np.asarray(s["xyz"], float) for k, s in enumerate(scene_specs)
            if k != mobile_idx and s.get("sym") is not None and len(s["sym"])]
    rest += [seg["ca"] for k, s in enumerate(scene_specs) if k != mobile_idx
             and s.get("cartoon") for seg in s["cartoon"]["segments"]]
    rest_c = np.vstack(rest).mean(0) if rest else mob_xyz.mean(0) - 1
    if along is None:
        along = mob_xyz.mean(0) - rest_c
    along = np.asarray(along, float)
    along /= (np.linalg.norm(along) + 1e-9)

    full = np.vstack(rest + [mob_xyz]) if rest else mob_xyz
    frames = []
    # camera frame fixed for the whole in-slide
    cam0 = _auto_view(full, view=view)

    def scene(dvec, spin_deg=0.0):
        pl = pv.Plotter(off_screen=True, window_size=size, lighting="none")
        pl.set_background(bg)
        for k, sp in enumerate(scene_specs):
            _draw_spec(pl, sp, offset=dvec if k == mobile_idx else None)
        _light_rig(pl)
        try:
            pl.enable_ssao(radius=2.0, bias=0.5, kernel_size=64, blur=True)
        except Exception:
            pass
        pl.camera_position = list(cam0)
        pl.camera.zoom(zoom)
        if spin_deg:
            pl.camera.Azimuth(spin_deg)
        try:
            pl.enable_anti_aliasing("msaa")
        except Exception:
            pass
        im = np.asarray(pl.screenshot(return_img=True))[..., :3]
        pl.close()
        return im

    for i in range(n_in):
        t = i / (n_in - 1)
        ease = t * t * (3 - 2 * t)               # smoothstep
        d = along * approach * (1 - ease)
        frames.append(scene(d))
    frames += [frames[-1]] * n_hold
    for i in range(n_spin):
        frames.append(scene(np.zeros(3), spin_deg=360.0 * i / n_spin))
    return _write_mp4(frames, out_mp4, fps)

"""
mcc — multi-coiled-coil sequence generator.

You normally only need these helpers (see the notebook):
    mcc.sequence(spec)      -> amino-acid sequence for one architecture
    mcc.describe(spec)      -> print a readable summary (modules, loops, contigs)
    mcc.draw_topology(spec) -> helical-wheel topology diagram (orientation + loops)
    mcc.contigs(spec)       -> MASTER contigs string (fixed helices; loops are the gaps)
    mcc.save_fasta(specs)   -> write all sequences to a FASTA file
    mcc.fold(...)           -> run local ColabFold on the FASTA (Linux + GPU)
    mcc.view_models(...)    -> browse the 5 folded models, pick one, save it

Spec grammar:  '6(2)4(2)6', '4(2)4', '-4(2)4', '3(1)3'
    numbers      = helices per module (2, 3, 4 or 6)
    (n)          = helices shared with the next module (1 or 2)
    +/- on a mod = handedness of module 1 only (+/none = CW, - = ACW); the rest propagate.

Advanced users: retune the chemistry in seq_lib / loop_lib below.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
import re, os, glob, json, shutil, subprocess

# ============================================================
# RULES  (edit these to retune the chemistry)
# ============================================================
REGISTRY = "bcdefga"        # residues stored & helices written N->C in this order (helix starts at b)
HEPTADS_PER_HELIX = 4

# heptad blocks per oligomeric state (= helices in the module): (acidic_block, basic_block)
seq_lib = {
    2: ("EELKQEI", "KKLKQEI"),
    3: ("EELKQEL", "KKLKQEL"),
    4: ("EEIAKQL", "KKIAEQL"),
    6: ("EEIAKLL", "KKIAELL"),
}

loop_lib = {
    "ae":   "GG",    # intra-module loop over an a/e-a/e interface (CW starts here)
    "dg":   "KKGGL", # intra-module loop over a d/g-d/g interface (ACW starts here)
    "link": "KGG",   # module-to-module linker
}

# offset maps {position_in_shared_helix : role_in_next_module}
OFFSET_TWO = {                       # two shared helices -> register shift +1 / -1
    "minus1": ({"e": "d", "b": "a", "f": "e"}, -1),
    "plus1":  ({"g": "a", "c": "d", "f": "g"}, +1),
}
TWO_SHARED_SELECT = {                # (handedness, parity) of preceding module -> map
    ("CW", "even"): "minus1", ("ACW", "even"): "plus1",
    ("CW", "odd"):  "plus1",  ("ACW", "odd"):  "minus1",
}
OFFSET_ONE = {                       # one shared helix -> +2 / -2 (sign deferred, stubbed)
    "plus2": ({"f": "d", "c": "a"}, +2),
    "minus2": ({"b": "d", "f": "a"}, -2),
}
ONE_SHARED_DEFAULT = "plus2"

# ============================================================
# DATA MODEL
# ============================================================
@dataclass
class Module:
    n_helices: int
    handedness: str = "CW"
    register_shift: int = 0
    offset_map: Dict[str, str] = field(default_factory=dict)

    @property
    def parity(self) -> str:
        return "even" if self.n_helices % 2 == 0 else "odd"


@dataclass
class Topology:
    modules: List[Module]
    shares: List[int]

    @classmethod
    def from_string(cls, spec: str) -> "Topology":
        return parse_topology(spec)


@dataclass
class Segment:
    kind: str          # "helix" | "loop"
    label: str
    seq: str
    note: str = ""


# ============================================================
# HELPERS + PARSER
# ============================================================
def _flip(h): return "ACW" if h == "CW" else "CW"
def _start_interface(h): return "ae" if h == "CW" else "dg"


def interface_types(mod: Module) -> List[str]:
    start = _start_interface(mod.handedness)
    other = "dg" if start == "ae" else "ae"
    return [start if i % 2 == 0 else other for i in range(mod.n_helices - 1)]


def parse_topology(spec) -> Topology:
    if not isinstance(spec, str):
        raise TypeError(f"expected a spec string like '6(2)4(2)6', got {type(spec).__name__}: {spec!r}")
    spec = spec.replace(" ", "")
    tokens = re.findall(r"[+-]?\d+[+-]?|\(\d+\)", spec)
    raw_mods, shares = [], []
    for t in tokens:
        if t.startswith("("):
            shares.append(int(t[1:-1]))
        else:
            n = int(re.search(r"\d+", t).group())
            sign = "-" if "-" in t else ("+" if "+" in t else "")
            raw_mods.append((n, sign))

    modules = [Module(n_helices=n) for n, _ in raw_mods]
    modules[0].handedness = "ACW" if raw_mods[0][1] == "-" else "CW"

    for i in range(len(modules) - 1):
        cur, nxt, sh = modules[i], modules[i + 1], shares[i]
        if sh == 2:
            omap, shift = OFFSET_TWO[TWO_SHARED_SELECT[(cur.handedness, cur.parity)]]
            # domain crossing adds one flip: flip on even, keep on odd
            nxt.handedness = _flip(cur.handedness) if cur.parity == "even" else cur.handedness
        elif sh == 1:
            omap, shift = OFFSET_ONE[ONE_SHARED_DEFAULT]
            nxt.handedness = cur.handedness            # stub: deferred +-2 sign
        else:
            raise ValueError(f"shared count must be 1 or 2, got {sh}")
        nxt.register_shift = shift
        nxt.offset_map = omap
    return Topology(modules=modules, shares=shares)


# ============================================================
# SEQUENCE GENERATION
# ============================================================
def _heptad_blocks(n_state, n_heptads=HEPTADS_PER_HELIX):
    acidic, basic = seq_lib[n_state]
    half = n_heptads // 2
    return [acidic] * half + [basic] * (n_heptads - half)


def _build_helix(n_state):
    return "".join(_heptad_blocks(n_state))


def _apply_offset(helix_seq, next_n_state, offset_map):
    """Shared helix: base = preceding module; offset positions carry the next module's
    g/a/d/e residues so both bundles are satisfied."""
    if not offset_map:
        return helix_seq
    acidic, _ = seq_lib[next_n_state]
    out = []
    for j in range(0, len(helix_seq), 7):
        hept = list(helix_seq[j:j + 7])
        for pos, role in offset_map.items():
            hept[REGISTRY.index(pos)] = acidic[REGISTRY.index(role)]
        out.append("".join(hept))
    return "".join(out)


def generate(spec) -> List[Segment]:
    """spec string (or Topology) -> ordered list of helix/loop Segments."""
    topo = spec if isinstance(spec, Topology) else Topology.from_string(spec)
    segments, M = [], topo.modules
    for mi, mod in enumerate(M):
        shared_prev = topo.shares[mi - 1] if mi > 0 else 0
        shared_next = topo.shares[mi] if mi < len(topo.shares) else 0
        itypes = interface_types(mod)
        for hi in range(mod.n_helices):
            if hi < shared_prev:
                continue
            if segments:
                if hi == shared_prev and mi > 0:
                    segments.append(Segment("loop", "link", loop_lib["link"]))
                else:
                    t = itypes[hi - 1]
                    segments.append(Segment("loop", t, loop_lib[t]))
            s = _build_helix(mod.n_helices)
            note = ""
            if shared_next and hi >= mod.n_helices - shared_next:
                nxt = M[mi + 1]
                s = _apply_offset(s, nxt.n_helices, nxt.offset_map)
                note = f"shared w/ M{mi + 2}"
            segments.append(Segment("helix", f"M{mi + 1}H{hi + 1}", s, note))
    return segments


def sequence(spec) -> str:
    """Full amino-acid sequence for an architecture."""
    return "".join(s.seq for s in generate(spec))


def annotate(spec):
    """[(Segment, start, end)] with 1-based residue ranges (chain A)."""
    out, pos = [], 1
    for s in generate(spec):
        end = pos + len(s.seq) - 1
        out.append((s, pos, end))
        pos = end + 1
    return out


def contigs(spec, chain="A") -> str:
    """MASTER contigs = fixed helix ranges; loops are the gaps between them."""
    return "/".join(f"{chain}{a}-{b}" for s, a, b in annotate(spec) if s.kind == "helix")


# ============================================================
# FRIENDLY OUTPUT + FILES
# ============================================================
def describe(spec) -> None:
    """Print a readable summary of one architecture."""
    topo = Topology.from_string(spec)
    segs = generate(topo)
    seq = "".join(s.seq for s in segs)
    n_hx = sum(1 for s in segs if s.kind == "helix")
    print(f"{spec}   ({len(seq)} aa, {n_hx} helices)")
    print("   modules :", ", ".join(f"M{i+1}={m.n_helices}{m.handedness}"
                                     for i, m in enumerate(topo.modules)))
    print("   loops   :", " ".join(s.label for s in segs if s.kind == "loop"))
    print("   contigs :", contigs(spec))
    print("   sequence:", seq)


def save_fasta(specs, path="mcc_sequences.fasta", width=60) -> str:
    """Write all architectures to a FASTA file (header = spec)."""
    with open(path, "w") as f:
        for spec in specs:
            s = sequence(spec)
            f.write(f">{spec}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + "\n")
    print(f"wrote {len(specs)} sequences -> {path}")
    return path


# ============================================================
# FOLDING (local ColabFold) + VIEWER
# ============================================================
def fold(fasta="mcc_sequences.fasta", results="results",
         msa_mode="single_sequence", extra=None) -> None:
    """Run a local colabfold_batch (needs a GPU + colabfold_batch on PATH).
    Equivalent to:  colabfold_batch --msa-mode single_sequence <fasta> <results>"""
    exe = shutil.which("colabfold_batch")
    if exe is None:
        raise RuntimeError("colabfold_batch not found on PATH. Install localcolabfold "
                           "and/or activate its environment, then restart the kernel.")
    cmd = [exe, "--msa-mode", msa_mode, fasta, results] + (extra or [])
    print("running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("done ->", results + "/")


def draw_topology(spec, save=None, add_labels=False, ax=None):
    """Topology diagram, drawn by the project's topology_utils.py (your convention).

    Modules are stacked as polygons (hexagon=6, square=4, triangle=3, pair=2); the chain
    walks each polygon and the last two helices of a module are its shared edge with the
    next. Circles are coloured by which terminus faces you (blue = N, red = C), with the
    chain's first/last helix ringed to mark the sequence N-/C-termini. Background is
    transparent when saved. Pass save="name.png"/".svg" to write a file.
    """
    import matplotlib.pyplot as plt
    import topology_utils as tu

    topo = Topology.from_string(spec)
    polygons = [m.n_helices for m in topo.modules]
    shared   = topo.shares
    windings = [m.handedness for m in topo.modules]    # CW/ACW per domain -> winding
    angles   = [[360.0 / n] * n for n in polygons]     # regular-polygon turn angles
    combined = tu.combine_domain_angles(angles, shared, windings=windings)
    x, y = tu.calc_cart_coords(combined)

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 8))
    tu.create_plot(x, y, ax=ax, add_labels=add_labels, outpath=save)


def _list_models(results):
    jobs = {}
    for p in glob.glob(os.path.join(results, "*_rank_*.pdb")):
        m = re.match(r"(.+?)_(unrelaxed|relaxed)_rank_(\d+)_", os.path.basename(p))
        if m:
            job, kind, rank = m.group(1), m.group(2), int(m.group(3))
            jobs.setdefault(job, {}).setdefault(rank, {})[kind] = p
    return jobs


def _mean_plddt(job, rank, results):
    f = glob.glob(os.path.join(results, f"{job}_scores_rank_{rank:03d}_*.json"))
    if not f:
        return None
    d = json.load(open(f[0]))
    return sum(d["plddt"]) / len(d["plddt"])


def view_models(results="results", save_dir="selected"):
    """Interactive: pick a sequence + model (1-5), view coloured by pLDDT, save your choice."""
    jobs = _list_models(results)
    if not jobs:
        print(f"No models found in '{results}/'. Run mcc.fold(...) first.")
        return
    py3Dmol = __import__("py3Dmol")
    from ipywidgets import Dropdown, Button, HBox, VBox, Output, HTML
    from IPython.display import display

    def path_for(job, rank):
        d = jobs[job][rank]
        return d.get("relaxed", d.get("unrelaxed"))

    job_dd = Dropdown(options=sorted(jobs), description="sequence")
    rank_dd = Dropdown(options=sorted(jobs[sorted(jobs)[0]]), description="rank")
    save_btn = Button(description="⭐ Save this model", button_style="success")
    status, view = HTML(), Output()

    def refresh(*_):
        with view:
            view.clear_output(wait=True)
            p = path_for(job_dd.value, rank_dd.value)
            pl = _mean_plddt(job_dd.value, rank_dd.value, results)
            print(f"{os.path.basename(p)}   mean pLDDT = {pl:.1f}" if pl else os.path.basename(p))
            v = py3Dmol.view(width=520, height=420)
            v.addModel(open(p).read(), "pdb")
            v.setStyle({"cartoon": {"colorscheme":
                       {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}})
            v.zoomTo(); v.show()

    def on_job(*_):
        rank_dd.options = sorted(jobs[job_dd.value]); refresh()

    def on_save(_):
        os.makedirs(save_dir, exist_ok=True)
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", job_dd.value)
        dst = os.path.join(save_dir, f"{safe}_rank{rank_dd.value}.pdb")
        shutil.copy(path_for(job_dd.value, rank_dd.value), dst)
        status.value = f"saved &rarr; <b>{dst}</b>"

    job_dd.observe(on_job, "value"); rank_dd.observe(refresh, "value")
    save_btn.on_click(on_save)
    display(VBox([HBox([job_dd, rank_dd, save_btn]), status, view]))
    refresh()

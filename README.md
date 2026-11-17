[Uploading README_mCC_generator_1_0.md…]()
# mCC generator 1.0 — User Guide

**mCC generator** is a single-file, offline tool for designing **multi-coiled-coil (mCC)**
proteins. You draw the topology — helices grouped into coiled-coil **domains**, joined into
one **chain** — and it generates the amino-acid **sequence** live, following coiled-coil
heptad rules. No install, no server, no internet.

- **File:** `mCC_generator_1_0.html`
- **To open:** double-click it (or drag it into Chrome, Edge, Safari or Firefox). It runs
  entirely in the browser; nothing is uploaded anywhere.
- **Folding is separate:** this tool only *designs* and writes the sequence. To fold a
  design with ColabFold on a GPU, use the companion server tool (see the folding-server
  README).

---

## Concepts (30 seconds)

- **Helix** — one α-helix, drawn as a circle.
- **Domain** (a coiled-coil bundle) — a ring of **2, 3, 4 or 6** helices (the sequence
  library also covers 5, 7, 8). The core-facing **a** and **d** positions of the heptad
  (`b c d e f g a`) pack the bundle together.
- **Shared helix** — a helix that belongs to **two or three** domains (shown with a gold
  dashed ring). This is how domains connect. Sharing **2 helices** = a shared *edge*;
  sharing **1 helix** = a shared *vertex* (a Y-junction).
- **Chain** — the order the helices are connected N→C by loops. The chain defines the
  final single-polypeptide sequence.

---

## The workflow

The controls are in the **left panel**, colour-coded by step.

### ① Build domains

Two ways, both starting from a selection:

1. Press **Select helices** and click helices on the canvas to highlight them (click again
   to un-highlight).
2. Then press **one** of:
   - **Add domain** — creates a **new** domain of the chosen **shape** (pick it from the
     shape dropdown) that is *connected to* the highlighted helices. Highlight **nothing**
     first to drop a standalone domain. Share **1 or 2** helices to connect it to existing
     structure (a helix in 2–3 domains becomes a shared helix).
   - **Group selected** — turns the highlighted **existing** helices into a domain (adds no
     new helices).

**Edit helices** (hand placement): **＋ Helix** (or double-click empty canvas), **Move**
(drag), **Delete** (click a helix, or right-click any helix). Remove a whole domain with
**✕ remove** in the *Domains* panel on the right. **Ctrl/⌘-Z** undoes.

### ② Chain order (N→C)

Press **Set N→C**, then click helices in order from the N-terminus to the C-terminus.
Arrows show the path; the first and last helix get an N/C ring. Click a helix again to
remove it from the chain.

### Read the sequence

The **Sequence** panel (right) updates live as you build — one line per helix, spaces
marking each heptad boundary, loops shown in orange. The status line under the canvas tells
you what's missing ("chain 3/4 set", etc.) and shows **"✓ chain complete"** when every
helix is in the chain. Use **Download FASTA** or **Copy**.

### Outputs

- **Download FASTA** — the sequence as a `.fasta` file (one unbroken chain).
- **Save PNG** (Tools) — the topology diagram you drew, as a transparent PNG.
- **Idealized** (Tools) — a clean, idealized schematic of the topology, which you can also
  save as PNG.

---

## What the colours and marks mean

- **Blue circle** = you're looking down that helix's **N-terminus**; **red** = its
  **C-terminus** (alternates along the chain).
- **Ringed circle** = the chain's **start (N)** and **end (C)**.
- **Gold dashed ring** = a **shared helix** (double ring = shared by 3 domains).
- **Dashed orange connector** = a **long loop**: the chain jumps between two helices that
  share no domain. It still generates (with a short `KGG` linker) but is flagged — check it
  makes geometric sense.

---

## Loops

Loops between helices are chosen automatically from the geometry: `GG` (a/e interface),
`KKGGL` (d/g interface), or `KGG` (a domain jump / after shared helices). You can override
any loop by typing into the **Loops** panel; press **auto** to revert.

---

## Shared helices and Y-junctions (notes)

- A canonical heptad helix has two core faces (a/e and d/g), so a helix can comfortably
  bridge **two** domains. The tool also allows a **three-domain** junction (one helix giving
  its a/d core to one domain and its a/e and d/g flanks to the other two). Junctions that
  over-subscribe a face (two domains competing for the same flank, or a domain forced onto
  the solvent side) are **flagged under the Sequence panel** — treat these as hypotheses to
  test by folding.

---

## Known limitation — the Idealized diagram

The **Idealized** schematic can only lay out designs where **every domain connects to the
rest through a shared 2-helix edge** (i.e. a linear or edge-glued chain of bundles).
**Branched** topologies (a hub with several spokes) and **single-shared-helix
(Y-junction)** designs can't be auto-laid-out, so for those the Idealized window shows a
short message instead of a picture.

This affects **only** the idealized schematic. The main canvas drawing, the live sequence,
the FASTA export, and folding all work normally for those designs — you just won't get the
tidy idealized picture.

---

## Supported bundle sizes

The shape dropdown offers **2, 3, 4, 6**. The sequence library also contains **5, 7, 8**, so
if you manually **Group selected** exactly 5, 7, or 8 helices into a domain it will still
generate. Any other size produces a clear message in the Sequence panel instead of a bad
sequence.

---

## Tips

- Start empty: pick a shape and press **Add domain** to drop your first bundle, then build
  outward by selecting its helices and adding more.
- **Esc** returns to Move mode; **Ctrl/⌘-Z** / **Shift-Ctrl/⌘-Z** undo / redo.
- Everything lives in the browser tab — there's no autosave, so **Download FASTA** (and
  Save PNG) before you close.

Questions about the chemistry — shapes, shared helices, a/d packing, loops — ask Kasia
(katarzyna.ozga@bio.ku.dk).

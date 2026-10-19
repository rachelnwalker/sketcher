# Template comparison gallery

Generate 2–5 distinct molecules per public ring-template reference and compare
fresh coordinates with `useRingTemplates=True` and `False` on a specified native
Sketcher/Core Suite runtime. The checked-in gallery was generated using 26-4.

## Generate locally

From the Sketcher repository:

```sh
python3 tools/template_gallery/run_local.py \
  --build-root /path/to/26-4/build \
  --references tools/template_gallery/references.smi \
  --sketcher-repo . \
  --runtime-label 26-4 \
  --output docs/template-gallery
python3 tools/template_gallery/validate.py docs/template-gallery
python3 tools/template_gallery/run_local.py --build-root /path/to/26-4/build --test
python3 -m http.server 8765 --directory docs
```

Open `http://localhost:8765/template-gallery/`. A current browser with
`DecompressionStream` support is needed for the small, compressed per-template
datasets. The page makes no requests to drawing services or external CDNs.

The launcher finds the build's own Python, RDKit package and native dependencies;
it does not install a different RDKit from PyPI. The generator can also be run
directly with an already configured build Python. For a quick development sample,
pass `--limit 8`; the manifest and page mark such output as incomplete, and the
deployment validator rejects it.

## Collection and provenance

`references.smi` is a snapshot of the public
[rdkit/molecular_templates](https://github.com/rdkit/molecular_templates) corpus.
Its upstream revision and checksum are in `references-source.json`. Its license is
in `REFERENCES-LICENSE`. Reference entries with different source layouts remain
separate even when they describe the same molecule.

The generator reads template CXSMARTS from the actual native Depictor library and
matches every reference against that compiled collection using atom count and
the ordered source coordinates. A missing reference is an error. The one extra
public cage template in Sketcher's patch is instantiated as an all-carbon
reference with its supplied coordinates, and is labeled as such. This gives 582
templates for the current 26-4 runtime, including all 581 reference entries.

The manifest records the measured runtime version and library SHA-256 separately
from the Sketcher commit, dependency pins and patch hashes. Checking in a new
recipe does not silently relabel old generated results as a new build. Regenerate
the dataset whenever the runtime changes.

## Molecule generation and comparison

Each template has its reference molecule plus up to four valid sidechain
variants: a short alkyl, a bulky/aryl group, a polar group, and a double
substitution. Attachment sites prioritize bridgeheads and distinct symmetry
classes. Variants are sanitized and deduplicated by isomeric canonical SMILES.
If a reference has no replaceable hydrogen, the generator may reduce one double
bond while preserving the template graph; this is labeled in the result. A
template that cannot yield at least two molecules remains visible as a coverage
limitation, and validation fails rather than silently omitting it.

Where possible, attachment sites retain a match to the source template SMARTS.
Cases for which this is impossible are labeled "Reference template ineligible";
another template in the full library may still match.

Both sides use the same molecule and atom order, remove existing conformers, force
the RDKit backend, and use deterministic settings. The only coordinate-generation
parameter that differs is `useRingTemplates`. The enabled side uses the runtime's
complete default library: it is not forced to use the reference template shown
in the header. RDKit may choose another matching template or use no template.

After generation, the disabled scaffold is aligned to the enabled scaffold with
rotation and translation only. Both SVGs use the same canvas bounds and scale.
Sidechain atoms are highlighted. RMSD uses the original atom correspondence;
no symmetry remapping, reflection, or scaling is used.

Close contacts, proper bond crossings, sidechain atoms inside scaffold ring
polygons, and specified-stereo MolBlock round trips are inspection aids. Crossings
can be intentional in cage drawings. The flags are not automatic quality grades.
Failures are retained in the manifest and gallery. `results.csv` contains all
comparison metrics without SVG data.

## Checks and publication

`validate.py` uses only the Python standard library and checks the complete
collection, unique molecule counts, both SVGs, canvas agreement, finite metrics,
CSV coverage and source provenance. Chemistry/geometry unit checks are in
`test/template_gallery/test_generation.py` and require the same RDKit runtime.

The existing GitHub Pages workflow validates and copies `docs/template-gallery`
to `/template-gallery/` alongside the WebAssembly Sketcher. It publishes the
checked-in 26-4 results and does not regenerate them with an unrelated CI runtime.

Gallery-only pushes use `deploy-gallery.yml` to reuse the app artifact from the
latest successful Pages app build. Mixed app/gallery changes use the full build.
If app artifacts have expired, run the full Pages workflow to rebuild them.

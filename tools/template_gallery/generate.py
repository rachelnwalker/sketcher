#!/usr/bin/env python3
"""Generate a public template/sidechain comparison using the active RDKit build.

Only a CXSMILES template corpus is accepted. No molecule uploads or remote calls
are made. Each output record contains both depictions of the same atom order.
"""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
from rdkit import Chem, Geometry, RDLogger, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

FRAGMENTS = {
    'methyl': 'C',
    'ethyl': 'CC',
    'isopropyl': 'C(C)C',
    'tert-butyl': 'C(C)(C)C',
    'phenyl': 'c1ccccc1',
    'hydroxy': 'O',
    'amino': 'N',
    'methoxy': 'OC',
    'amide chain': 'CC(=O)N',
}
WIDTH, HEIGHT = 500, 340


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def revision(path):
    result = subprocess.run(
        ['git', '-C', str(path), 'rev-parse', 'HEAD'],
        capture_output=True,
        text=True,
        check=False)
    return result.stdout.strip() if result.returncode == 0 else None


def coordinate_signature(mol):
    return (mol.GetNumAtoms(),
            tuple(
                round(float(v), 3)
                for row in mol.GetConformer().GetPositions()
                for v in row[:2]))


def find_library():
    import rdkit
    lib_dir = Path(rdkit.__file__).resolve().parent.parents[2]
    choices = list(lib_dir.glob('libRDKitDepictor*.dylib'))
    choices += list(lib_dir.glob('libRDKitDepictor*.so*'))
    if not choices:
        raise ValueError(
            'Pass --template-library pointing at this build\'s Depictor library'
        )
    return sorted(choices)[0].resolve()


def load_templates(references, library):
    # The literal CXSMARTS are read from the actual native library, rather than
    # assuming the latest upstream corpus is the one compiled into the build.
    queries = []
    for literal in library.read_bytes().split(b'\0'):
        if literal.startswith(b'[!#200') and b' |(' in literal:
            value = literal.decode('utf-8')
            query = Chem.MolFromSmarts(value)
            if query is not None and query.GetNumConformers():
                queries.append((value, query))
    if not queries:
        raise ValueError(
            'No built-in CXSMARTS found in the supplied Depictor library')
    signatures = {coordinate_signature(m): (s, m) for s, m in queries}
    records, seen = [], set()
    for line_number, line in enumerate(references.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        mol = Chem.MolFromSmiles(line)
        if mol is None or not mol.GetNumConformers():
            raise ValueError(
                f'Invalid reference CXSMILES on line {line_number}')
        key = coordinate_signature(mol)
        if key not in signatures:
            raise ValueError(
                f'Reference line {line_number} is absent from the runtime library'
            )
        seen.add(key)
        records.append({
            'mol': mol,
            'cxsmiles': line,
            'line': line_number,
            'referenceKind': 'public corpus',
            'query': signatures[key][0]
        })
    # The Sketcher patch includes one extra, public, all-carbon cage topology.
    # Instantiate its query graph as carbon and keep its supplied coordinates.
    # Never substitute an unrelated private molecule as a missing reference.
    for key, (query_text, query) in signatures.items():
        if key in seen:
            continue
        concrete = re.sub(r'\[[^\]]*!#200[^\]]*\]', 'C',
                          query_text).replace('~', '-')
        mol = Chem.MolFromSmiles(concrete)
        if mol is None or not mol.GetNumConformers():
            raise ValueError(
                'A runtime-only template needs a chemically valid reference molecule'
            )
        records.append({
            'mol': mol,
            'cxsmiles': concrete,
            'line': None,
            'referenceKind': 'all-carbon reference from public runtime template',
            'query': query_text
        })
    # Duplicate template literals may be pooled by a C++ compiler. Do not claim
    # an upstream record absent just because the binary stores a string once.
    return records, len(signatures)


def canonical(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def reference_matches(mol, query, core_count):
    target = Chem.RWMol(mol)
    for atom in target.GetAtoms():
        if atom.GetIdx() >= core_count:
            atom.SetAtomicNum(200)
    return target.GetMol().HasSubstructMatch(query)


def sites(mol):
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    atom_rings = mol.GetRingInfo().AtomRings()
    choices = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() > 1 and a.GetTotalNumHs() > 0
    ]
    # Cover bridgeheads before less crowded sites, then use symmetry classes
    # so four variants don't repeatedly exercise the same environment.
    choices.sort(key=lambda i: (-sum(i in r for r in atom_rings), -mol.
                                GetAtomWithIdx(i).GetDegree(), ranks[i], i))
    unique, seen = [], set()
    for i in choices:
        if ranks[i] not in seen:
            unique.append(i)
            seen.add(ranks[i])
    return unique + [i for i in choices if i not in unique]


def attach(mol, idx, fragment):
    group = Chem.MolFromSmiles(fragment)
    offset = mol.GetNumAtoms()
    result = Chem.RWMol(Chem.CombineMols(mol, group))
    atom = result.GetAtomWithIdx(idx)
    if atom.GetNumExplicitHs():
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
    result.AddBond(idx, offset, Chem.BondType.SINGLE)
    result.RemoveAllConformers()
    Chem.SanitizeMol(result)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    return Chem.Mol(result), list(range(offset, result.GetNumAtoms()))


def reduced_scaffold(mol):
    # A template's SMARTS matches any bond order. For a fully substituted
    # scaffold with no hydrogen, try reducing one bond to expose a valid site.
    rw = Chem.RWMol(mol)
    Chem.Kekulize(rw, clearAromaticFlags=True)
    for bond in rw.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        candidate = Chem.RWMol(rw)
        changed = candidate.GetBondWithIdx(bond.GetIdx())
        changed.SetBondType(Chem.BondType.SINGLE)
        changed.SetStereo(Chem.BondStereo.STEREONONE)
        for idx in (changed.GetBeginAtomIdx(), changed.GetEndAtomIdx()):
            candidate.GetAtomWithIdx(idx).SetNoImplicit(False)
        try:
            Chem.SanitizeMol(candidate)
            if sites(candidate):
                return Chem.Mol(candidate)
        except (ValueError, RuntimeError):
            continue
    return None


def make_variants(base, template_number, maximum, reference_query=None):
    parent = Chem.Mol(base)
    parent.RemoveAllConformers()
    output = [{
        'mol': parent,
        'label': 'Reference molecule',
        'substituent': 'none',
        'sites': [],
        'sidechainAtoms': [],
        'scaffoldChange': None,
        'referenceEligible': reference_matches(parent, reference_query,
                                               base.GetNumAtoms())
                             if reference_query is not None else None
    }]
    work, scaffold_change = parent, None
    available = sites(work)
    if not available:
        work = reduced_scaffold(parent)
        if work is None:
            return output, 'No chemically valid attachment site, including after one-bond reduction'
        scaffold_change = 'One ring double bond reduced to create an attachment site'
        available = sites(work)
    options = [
        ['methyl', 'ethyl'],
        ['isopropyl', 'tert-butyl', 'phenyl'],
        ['hydroxy', 'amino', 'methoxy', 'amide chain'],
        ['double substitution', 'phenyl', 'ethyl'],
    ]
    seen = {canonical(parent)}
    for slot, group in enumerate(options[:maximum - 1]):
        rotation = (template_number + slot) % len(group) if slot < 3 else 0
        names = group[rotation:] + group[:rotation]
        found, fallback = None, None
        for name in names:
            order = available[slot %
                              len(available):] + available[:slot %
                                                           len(available)]
            for idx in order:
                try:
                    if name == 'double substitution':
                        candidate, added = attach(work, idx, 'C')
                        others = [
                            i for i in sites(candidate)
                            if i < base.GetNumAtoms() and i != idx
                        ]
                        if not others and candidate.GetAtomWithIdx(
                                idx).GetTotalNumHs():
                            others = [idx]
                        if not others:
                            continue
                        candidate, more = attach(candidate, others[0], 'CC')
                        added += more
                        attachments = [idx, others[0]]
                    else:
                        candidate, added = attach(work, idx, FRAGMENTS[name])
                        attachments = [idx]
                    identity = canonical(candidate)
                    if identity in seen:
                        continue
                    eligible = reference_matches(
                        candidate, reference_query, base.GetNumAtoms(
                        )) if reference_query is not None else None
                    result = {
                        'mol': candidate,
                        'label': name.capitalize(),
                        'substituent': name,
                        'sites': attachments,
                        'sidechainAtoms': added,
                        'scaffoldChange': scaffold_change,
                        'referenceEligible': eligible
                    }
                    if eligible is False:
                        if fallback is None:
                            fallback = result
                        continue
                    found = result
                    break
                except (ValueError, RuntimeError):
                    continue
            if found:
                break
        found = found or fallback
        if found:
            seen.add(canonical(found['mol']))
            output.append(found)
    return output, None if len(
        output) >= 2 else 'No distinct valid sidechain variant'


def align_pair(enabled, disabled, core_count):
    target = np.array(enabled.GetConformer().GetPositions())[:, :2]
    points = np.array(disabled.GetConformer().GetPositions())[:, :2]
    source_center = points[:core_count].mean(axis=0)
    target_center = target[:core_count].mean(axis=0)
    u, _, vt = np.linalg.svd((points[:core_count] - source_center).T
                             @ (target[:core_count] - target_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    aligned = (points - source_center) @ rotation + target_center
    conf = disabled.GetConformer()
    for idx, (x, y) in enumerate(aligned):
        conf.SetAtomPosition(idx, Geometry.Point3D(float(x), float(y), 0))
    delta = np.sum((target - aligned)**2, axis=1)
    return float(np.sqrt(delta.mean())), float(
        np.sqrt(delta[:core_count].mean()))


def crossing(a, b, c, d):

    def cross(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

    # Proper intersections only: shared endpoints and collinear segments are
    # not counted as bond crossings.
    return cross(a, b, c) * cross(a, b, d) < -1e-9 and cross(c, d, a) * cross(
        c, d, b) < -1e-9


def inside(point, polygon):
    x, y = point
    result = False
    for p, q in zip(polygon, np.roll(polygon, -1, axis=0)):
        if (p[1] > y) != (
                q[1]
                > y) and x < (q[0] - p[0]) * (y - p[1]) / (q[1] - p[1]) + p[0]:
            result = not result
    return result


def geometry_metrics(mol, sidechain_atoms, core_count):
    points = np.array(mol.GetConformer().GetPositions())[:, :2]
    edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    bonded = {tuple(sorted(e)) for e in edges}
    lengths = [np.linalg.norm(points[a] - points[b]) for a, b in edges]
    scale = float(np.median(lengths)) if lengths else 1.5
    collisions = sum(
        np.linalg.norm(points[i] - points[j]) < 0.5 * scale
        for i in range(len(points))
        for j in range(i + 1, len(points))
        if (i, j) not in bonded)
    crossings = sum(
        crossing(points[a], points[b], points[c], points[d])
        for k, (a, b) in enumerate(edges)
        for c, d in edges[k + 1:]
        if len({a, b, c, d}) == 4)
    rings = [
        r for r in mol.GetRingInfo().AtomRings()
        if all(i < core_count for i in r)
    ]
    inward = sum(
        any(inside(points[i], points[list(r)])
            for r in rings)
        for i in sidechain_atoms)
    specified = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                    for a in mol.GetAtoms())
    specified |= any(b.GetStereo() not in (Chem.BondStereo.STEREONONE,
                                           Chem.BondStereo.STEREOANY)
                     for b in mol.GetBonds())
    roundtrip = None
    if specified:
        try:
            reread = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol),
                                          removeHs=False)
            roundtrip = reread is not None and canonical(reread) == canonical(
                mol)
        except (ValueError, RuntimeError):
            roundtrip = False
    return {
        'closeContacts': int(collisions),
        'bondCrossings': int(crossings),
        'sidechainAtomsInsideRings': int(inward),
        'stereoRoundTrip': roundtrip
    }


def draw_svg(mol, sidechain_atoms, bounds=None, width=WIDTH, height=HEIGHT):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height, -1, -1, True)
    opts = drawer.drawOptions()
    opts.prepareMolsBeforeDrawing = False
    opts.padding = 0.09
    opts.bondLineWidth = 1.7
    opts.minFontSize = 10
    opts.maxFontSize = 18
    opts.clearBackground = True
    opts.useBWAtomPalette()
    if bounds is not None:
        low, high = bounds
        drawer.SetScale(width, height, Geometry.Point2D(*map(float, low)),
                        Geometry.Point2D(*map(float, high)))
    prepared = rdMolDraw2D.PrepareMolForDrawing(mol,
                                                kekulize=True,
                                                addChiralHs=False,
                                                wedgeBonds=True,
                                                forceCoords=False)
    atoms = list(sidechain_atoms)
    color = {i: (1.0, 0.78, 0.46) for i in atoms}
    drawer.DrawMolecule(prepared,
                        highlightAtoms=atoms,
                        highlightAtomColors=color)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg[svg.index('<svg'):]
    # Fixed canvas coordinates are retained; the viewer adds responsive CSS.
    return svg


def depict_case(variant, core_count):
    enabled, disabled = Chem.Mol(variant['mol']), Chem.Mol(variant['mol'])
    for mol, templates in ((enabled, True), (disabled, False)):
        mol.RemoveAllConformers()
        rdDepictor.Compute2DCoords(mol,
                                   canonOrient=True,
                                   clearConfs=True,
                                   forceRDKit=True,
                                   useRingTemplates=templates,
                                   nFlipsPerSample=0,
                                   nSample=0,
                                   sampleSeed=0)
    rmsd, scaffold_rmsd = align_pair(enabled, disabled, core_count)
    points = np.vstack([
        np.array(m.GetConformer().GetPositions())[:, :2]
        for m in (enabled, disabled)
    ])
    low, high = points.min(axis=0) - 0.8, points.max(axis=0) + 0.8
    report = {k: v for k, v in variant.items() if k != 'mol'}
    report['sites'] = [i + 1 for i in variant['sites']]
    report.update(smiles=canonical(variant['mol']),
                  rmsd=round(rmsd, 5),
                  scaffoldRmsd=round(scaffold_rmsd, 5),
                  atomCount=enabled.GetNumAtoms())
    for label, mol in (('enabled', enabled), ('disabled', disabled)):
        report[label] = geometry_metrics(mol, variant['sidechainAtoms'],
                                         core_count)
        report[label]['svg'] = draw_svg(mol, variant['sidechainAtoms'],
                                        (low, high))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sketcher-repo', type=Path, required=True)
    parser.add_argument('--runtime-label', default='26-4')
    parser.add_argument('--template-library', type=Path)
    parser.add_argument('--max-molecules',
                        type=int,
                        choices=range(2, 6),
                        default=5)
    parser.add_argument(
        '--limit',
        type=int,
        help='Development sample only; manifest marks it incomplete')
    args = parser.parse_args()
    if args.references.suffix.lower() != '.smi':
        parser.error('--references must be a public .smi template corpus')
    library = args.template_library or find_library()
    RDLogger.DisableLog('rdApp.warning')
    RDLogger.DisableLog('rdApp.error')
    rdDepictor.LoadDefaultRingSystemTemplates()
    templates, runtime_count = load_templates(args.references, library)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = output / 'data'
    data.mkdir(exist_ok=True)
    summaries, errors, coverage, counts, csv_rows = [], [], [], Counter(), []
    selected = templates[:args.limit] if args.limit else templates
    for idx, record in enumerate(selected, 1):
        tid = f'T{idx:04d}'
        base = record['mol']
        reference_query = Chem.MolFromSmarts(record['query'])
        variants, limitation = make_variants(base, idx, args.max_molecules,
                                             reference_query)
        row = {
            'id': tid,
            'number': idx,
            'referenceSmiles': canonical(base),
            'referenceCXSmiles': record['cxsmiles'],
            'referenceKind': record['referenceKind'],
            'referenceLine': record['line'],
            'templateSmarts': record['query'],
            'atomCount': base.GetNumAtoms(),
            'cases': [],
            'limitation': limitation,
            'sourceSvg': draw_svg(base, [], width=240, height=150)
        }
        for variant in variants:
            try:
                row['cases'].append(depict_case(variant, base.GetNumAtoms()))
            except (ValueError, RuntimeError) as exc:
                failure = {
                    'label': variant['label'],
                    'smiles': canonical(variant['mol']),
                    'error': str(exc),
                    'substituent': variant['substituent']
                }
                row['cases'].append(failure)
                errors.append({'template': tid, **failure})
        usable = [c for c in row['cases'] if 'error' not in c]
        if len(variants) < 2:
            coverage.append({'template': tid, 'reason': limitation})
        flags = sum(
            any(c[mode][k]
                for k in ('closeContacts', 'bondCrossings',
                          'sidechainAtomsInsideRings')) or
            c[mode]['stereoRoundTrip'] is False
            for c in usable
            for mode in ('enabled', 'disabled'))
        summary = {
            'id': tid,
            'number': idx,
            'atomCount': base.GetNumAtoms(),
            'caseCount': len(row['cases']),
            'maxRmsd': max((c['rmsd'] for c in usable), default=0),
            'flags': flags,
            'errors': len(row['cases']) - len(usable),
            'limitation': limitation,
            'referenceSmiles': row['referenceSmiles'],
            'substituents': sorted({c['substituent'] for c in row['cases']}),
            'data': f'data/{tid}.json.gz'
        }
        summaries.append(summary)
        for case in row['cases']:
            flat = {
                'template': tid,
                'variant': case['label'],
                'smiles': case['smiles'],
                'reference_template_eligible': case.get('referenceEligible'),
                'rmsd': case.get('rmsd'),
                'scaffold_rmsd': case.get('scaffoldRmsd'),
                'error': case.get('error', '')
            }
            for mode in ('enabled', 'disabled'):
                for key in ('closeContacts', 'bondCrossings',
                            'sidechainAtomsInsideRings', 'stereoRoundTrip'):
                    flat[f'{mode}_{key}'] = case.get(mode, {}).get(key)
            csv_rows.append(flat)
        payload = json.dumps(row, separators=(',', ':'),
                             allow_nan=False).encode()
        (data / f'{tid}.json.gz').write_bytes(gzip.compress(payload, mtime=0))
        summary['dataSha256'] = sha256(data / f'{tid}.json.gz')
        counts['molecules'] += len(row['cases'])
        counts['successfulPairs'] += len(usable)
        counts['changedPairs'] += sum(c['rmsd'] > 0.05 for c in usable)
        if idx % 25 == 0 or idx == len(selected):
            print(
                f'{idx}/{len(selected)} templates; {counts["molecules"]} molecules',
                flush=True)
    patches = {
        p.name: sha256(p)
        for p in sorted((args.sketcher_repo / 'external').glob('rdkit*.patch'))
    }
    reference_revision = revision(args.references.parent)
    reference_metadata = args.references.with_name('references-source.json')
    if reference_metadata.exists():
        source = json.loads(reference_metadata.read_text())
        if source['sha256'] != sha256(args.references):
            raise ValueError(
                'Vendored reference corpus differs from its source metadata')
        reference_revision = source['commit']
    manifest = {
        'schemaVersion': 1,
        'runtimeLabel': args.runtime_label,
        'rdkitVersion': rdBase.rdkitVersion,
        'sketcherCommit': revision(args.sketcher_repo),
        'generatorSha256': sha256(Path(__file__)),
        'library': {
            'name': library.name,
            'sha256': sha256(library)
        },
        'dependencyPins': json.loads(
            (args.sketcher_repo / 'external/versions.json').read_text()),
        'recipePatchHashes': patches,
        'referenceRepository': 'https://github.com/rdkit/molecular_templates',
        'referenceCommit': reference_revision,
        'referenceSha256': sha256(args.references),
        'runtimeTemplateCount': runtime_count,
        'templateCount': len(selected),
        'complete': args.limit is None,
        'generatedAt': datetime.now(timezone.utc).isoformat(),
        'settings': {
            'useRingTemplates': [True, False],
            'forceRDKit': True,
            'canonOrient': True,
            'clearConfs': True,
            'sampleSeed': 0,
            'nSample': 0,
            'nFlipsPerSample': 0,
            'alignment': 'rotation and translation on the scaffold; no reflection or scaling',
            'maxMoleculesPerTemplate': args.max_molecules,
            'closeContactCutoff': '0.5 times the median bond length, nonbonded atom pairs'
        },
        'totals': dict(counts),
        'failures': errors,
        'coverageLimitations': coverage,
        'templates': summaries
    }
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    if csv_rows:
        with (output / 'results.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
            writer.writeheader()
            writer.writerows(csv_rows)
    assets = Path(__file__).resolve().parents[2] / 'docs/template-gallery'
    for name in ('index.html', 'gallery.js', 'gallery.css'):
        source = assets / name
        if source.exists() and source.resolve() != (output / name).resolve():
            shutil.copyfile(source, output / name)
    print(json.dumps({
        'templates': len(selected),
        **counts, 'coverageLimitations': len(coverage),
        'failures': len(errors)
    }),
          flush=True)
    return 0 if not coverage and not errors else 2


if __name__ == '__main__':
    sys.exit(main())

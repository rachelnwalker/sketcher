#!/usr/bin/env python3
"""Validate a generated gallery using only Python's standard library."""

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def validate(root, allow_sample=False):
    manifest = json.loads((root / 'manifest.json').read_text())
    assert manifest['schemaVersion'] == 1, 'Unknown dataset schema'
    assert allow_sample or manifest['complete'], 'Incomplete development sample'
    assert len(manifest['templates']) == manifest['templateCount']
    assert allow_sample or manifest['templateCount'] >= manifest[
        'runtimeTemplateCount']
    assert manifest['rdkitVersion'] and manifest['runtimeLabel']
    assert len(manifest['library']['sha256']) == 64
    assert manifest['referenceCommit'], 'Missing reference source revision'
    assert not manifest['failures'], 'The gallery has depiction failures'
    assert not manifest[
        'coverageLimitations'], 'Some templates lack valid sidechain coverage'
    settings = manifest['settings']
    assert settings['useRingTemplates'] == [True, False]
    assert settings['forceRDKit'] and settings['clearConfs']
    assert settings['sampleSeed'] == 0 and settings['nSample'] == 0
    pairs, changed = 0, 0
    ids = set()
    for item in manifest['templates']:
        assert item['id'] not in ids, 'Duplicate template ID'
        ids.add(item['id'])
        path = root / item['data']
        assert path.resolve().is_relative_to(root.resolve()), 'Unsafe data path'
        assert hashlib.sha256(path.read_bytes()).hexdigest(
        ) == item['dataSha256'], 'Dataset checksum mismatch'
        row = json.loads(gzip.decompress(path.read_bytes()))
        assert row['id'] == item['id']
        assert 2 <= len(
            row['cases']) <= 5, f'{row["id"]}: molecule count outside 2–5'
        assert len(row['cases']) == item['caseCount']
        identities = set()
        for case in row['cases']:
            assert 'error' not in case, f'{row["id"]}: failed depiction'
            assert case[
                'smiles'] not in identities, f'{row["id"]}: duplicate molecule'
            identities.add(case['smiles'])
            assert math.isfinite(case['rmsd']) and case['rmsd'] >= 0
            assert isinstance(case['referenceEligible'], bool)
            assert math.isfinite(case['scaffoldRmsd'])
            sizes = []
            for mode in ('enabled', 'disabled'):
                value = case[mode]
                drawing = ET.fromstring(value['svg'])
                assert drawing.tag.endswith('svg')
                assert not any(n.tag.endswith('script') for n in drawing.iter())
                sizes.append((drawing.attrib['width'], drawing.attrib['height'],
                              drawing.attrib.get('viewBox')))
                for key in ('closeContacts', 'bondCrossings',
                            'sidechainAtomsInsideRings'):
                    assert isinstance(value[key], int) and value[key] >= 0
            assert sizes[0] == sizes[1], 'Mismatched comparison canvases'
            pairs += 1
            changed += case['rmsd'] > 0.05
        assert row['cases'][0]['substituent'] == 'none', 'Missing baseline'
    assert pairs == manifest['totals']['molecules'] == manifest['totals'][
        'successfulPairs']
    assert changed == manifest['totals']['changedPairs']
    with (root / 'results.csv').open(newline='') as stream:
        assert sum(1 for _ in csv.DictReader(stream)) == pairs
    assert {p.stem.split('.')[0] for p in (root / 'data').glob('*.json.gz')
           } == ids
    for asset in ('index.html', 'gallery.js', 'gallery.css'):
        assert (root / asset).is_file(), f'Missing {asset}'
    return manifest['templateCount'], pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--allow-sample', action='store_true')
    args = parser.parse_args()
    templates, pairs = validate(args.directory, args.allow_sample)
    print(
        f'Validated {templates} templates, {pairs} unique molecules, {2 * pairs} SVG depictions'
    )


if __name__ == '__main__':
    main()

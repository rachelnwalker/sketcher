#!/usr/bin/env python3
"""Run the gallery generator with a Core Suite build's own Python and RDKit."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-root', type=Path, required=True)
    parser.add_argument('--test',
                        action='store_true',
                        help='Run the chemistry and geometry unit checks')
    args, extra = parser.parse_known_args()
    root = args.build_root.resolve()
    candidates = sorted(
        (root / 'schrodinger_buildenv_packages/.pixi/envs').glob('*'))
    environments = [
        p for p in candidates
        if list(p.glob('rdkit-*/lib/python*/site-packages/rdkit'))
    ]
    if len(environments) != 1:
        parser.error(
            'Expected one package environment containing RDKit; found ' +
            str(environments))
    packages = environments[0]
    rdkit = sorted(packages.glob('rdkit-*/lib/python*/site-packages'))
    python = sorted(packages.glob('python-*/bin/python3.*'))
    python = [
        p for p in python
        if p.name.count('.') == 1 and p.name.split('.')[-1].isdigit()
    ]
    if len(rdkit) != 1 or not python:
        parser.error('Could not identify the build Python and RDKit package')
    py_version = rdkit[0].parent.name
    sites = [rdkit[0], root / 'internal/lib' / py_version / 'site-packages']
    env = os.environ.copy()
    env['SCHRODINGER'] = str(root)
    mmshare = sorted(root.glob('mmshare-v*/bin/*'))
    if mmshare:
        env['MMSHARE_EXEC'] = str(mmshare[-1])
    env['PYTHONPATH'] = os.pathsep.join(str(p) for p in sites if p.is_dir())
    libraries = [packages / 'lib', *sorted(packages.glob('*/lib'))]
    env['DYLD_LIBRARY_PATH'] = os.pathsep.join(map(str, libraries))
    env['LD_LIBRARY_PATH'] = os.pathsep.join(map(str, libraries))
    script = (Path(__file__).resolve().parents[2] /
              'test/template_gallery/test_generation.py'
              if args.test else Path(__file__).with_name('generate.py'))
    return subprocess.call([str(python[-1]), str(script), *extra], env=env)


if __name__ == '__main__':
    sys.exit(main())

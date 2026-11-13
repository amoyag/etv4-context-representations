#!/usr/bin/env python3
"""Hash downloaded assets without deserialising weights or executing models."""
import hashlib
import json
from pathlib import Path
import platform

ROOT = Path(__file__).resolve().parents[1]
REV = 'fe4084f356480fd048a23097fcc133f1e5ba019d'
FILES = [
    ('ProFam-1', 'models/profam-1/checkpoints/last.ckpt',
     f'https://huggingface.co/judewells/ProFam-1/resolve/{REV}/checkpoints/last.ckpt',
     '8ec54d80323f188d3454fb4457f63f12158fe55552b919554806b12730c30e08', 1511154085),
    ('ProFam-1 config', 'models/profam-1/.hydra/config.yaml',
     f'https://huggingface.co/judewells/ProFam-1/resolve/{REV}/.hydra/config.yaml', None, 13415),
    ('ESM-2 650M', 'models/esm2/esm2_t33_650M_UR50D.pt',
     'https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt', None, 2604537549),
    ('ESM-2 contact regression', 'models/esm2/esm2_t33_650M_UR50D-contact-regression.pt',
     'https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt', None, 3687),
]

def main():
    records = []
    for model, relative, url, expected_hash, expected_size in FILES:
        path = ROOT / relative
        if not path.is_file():
            raise SystemExit(f'Missing asset: {relative}')
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        size = path.stat().st_size
        if not size or (expected_size is not None and size != expected_size):
            raise SystemExit(f'Invalid size: {relative}')
        if expected_hash and digest != expected_hash:
            raise SystemExit(f'Hash mismatch: {relative}')
        records.append(dict(model=model, path=relative, bytes=size, sha256=digest,
                            source=url, upstream_sha256=expected_hash,
                            upstream_revision=REV if model.startswith('ProFam') else None,
                            identity_check='upstream_sha256_match' if expected_hash else
                            'local_sha256_recorded_from_official_download'))
    report = dict(scope='asset_download_and_hash_only', host=platform.node(),
                  recorded_date='2026-11-10', files=records, passed=True,
                  weights_loaded=False, inference_executed=False, E1_qualified=False,
                  limitation='ESM hashes record downloaded identity; no independent upstream SHA-256 was available.')
    (ROOT / 'environment/model-assets-macmini.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()

"""Offline allowlist and credential/path-pattern check. Never print matched secrets."""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
ALLOWED = {'.R', '.md', '.json', '.py'}
PATTERNS = {
    'github_token': r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,})\b',
    'jwt_literal': r'\beyJ[A-Za-z0-9_-]{15,}\.[A-Za-z0-9_-]{15,}\.[A-Za-z0-9_-]{10,}\b',
    'private_key': r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
    'personal_windows_path': r'[A-Za-z]:[/\\]+Users[/\\]+',
}
errors = []
files = []
for path in sorted(ROOT.rglob('*')):
    rel = path.relative_to(ROOT)
    if '.git' in rel.parts or '__pycache__' in rel.parts or not path.is_file():
        continue
    files.append(path)
    if path.name not in {'.gitignore', '.gitattributes'} and path.suffix not in ALLOWED:
        errors.append(f'Unapproved file type: {rel}')
        continue
    text = path.read_text(encoding='utf-8')
    if '\ufffd' in text or '\x00' in text:
        errors.append(f'Encoding failure: {rel}')
    for label, pattern in PATTERNS.items():
        if re.search(pattern, text):
            errors.append(f'{label}: {rel}')
for item in json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text(encoding='utf-8')):
    digest = hashlib.sha256((ROOT / item['file']).read_bytes()).hexdigest()
    if digest != item['release_sha256']:
        errors.append(f'Checksum mismatch: {item["file"]}')
if errors:
    raise SystemExit('\n'.join(errors))
print(f'PASS: {len(files)} code/documentation files; source checksums match; no flagged credential or personal-path patterns.')
print('This automated scan complements, but does not replace, manual inspection.')

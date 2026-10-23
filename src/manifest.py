"""Write data_manifest.csv: every raw input with its byte length and a fingerprint.
Files up to 32 MiB get a full sha256; larger ones a headtail256 (sha256 over the first
and last mebibyte plus the exact length), the convention used across these projects."""
import csv
import hashlib

from common import RAW, ROOT

MIB = 1 << 20


def fingerprint(path):
    size = path.stat().st_size
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        if size <= 32 * MIB:
            h.update(fh.read())
            return size, "sha256:" + h.hexdigest()
        h.update(fh.read(MIB))
        fh.seek(size - MIB)
        h.update(fh.read(MIB))
        h.update(str(size).encode())
        return size, "headtail256:" + h.hexdigest()


rows = []
for p in sorted(RAW.rglob("*")):
    if p.is_file():
        size, fp = fingerprint(p)
        rows.append({"path": p.relative_to(ROOT).as_posix(), "bytes": size, "fingerprint": fp})
with open(ROOT / "data_manifest.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["path", "bytes", "fingerprint"])
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} files, {sum(r['bytes'] for r in rows) / 1e6:.1f} MB")

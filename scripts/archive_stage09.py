import pathlib
import zipfile
import hashlib
import csv

root = pathlib.Path(__file__).resolve().parents[1]
archive = root / "results/TNBC_ResponseResiduals_FinalAnalysis.zip"
folders = ["config", "metadata", "R", "scripts", "results", "tables", "figures", "logs"]
files = [root / "renv.lock"]
for folder in folders:
    files.extend(p for p in (root / folder).rglob("*") if p.is_file())
files = sorted(set(p for p in files if p != archive and p.suffix.lower() != ".rds"
                   and not p.name.endswith(".part")))
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for file in files:
        z.write(file, file.relative_to(root).as_posix())
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None, "Archive CRC validation failed"
    assert "metadata/data_inventory.tsv" in z.namelist()
    assert "renv.lock" in z.namelist()
    assert not any(n.startswith(("data/raw/", "models/")) or n.lower().endswith(".rds") for n in z.namelist())
print(f"Archive verified: {len(files)} files; {archive.stat().st_size} bytes")

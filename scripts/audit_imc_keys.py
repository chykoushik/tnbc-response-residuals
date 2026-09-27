import csv
import platform
from pathlib import Path

import h5py

root = Path("D:/TNBC-ResponseResiduals")
with (root / "metadata/neotrip_mapping_audit.tsv").open(newline="") as stream:
    identity = list(csv.DictReader(stream, delimiter="\t"))
keys = set()
for row in identity:
    for field in ("gsm", "geo_library_name", "geo_patient_id"):
        if row[field]:
            keys.add(row[field])
rows = []
with h5py.File(root / "external/neotrip_2026_paper/nt_preprocessed.h5ad", "r") as data:
    observations = data["obs"]
    for field in ("Patient_ID", "Tissue_ID", "fov", "_index"):
        node = observations[field]
        values = node["categories"] if isinstance(node, h5py.Group) else node
        if values.dtype.kind in "OSU":
            values = set(values.asstr()[:])
        else:
            values = set(str(value) for value in values[:])
        matches = keys.intersection(values)
        rows.append((field, len(values), len(matches), "Exact string comparison only; no prefix conversion, numeric matching or inferred crosswalk"))
    schema = [(name, "not_read_outcome_or_treatment" if name in ("Arm", "pCR") else "schema_only")
              for name in observations]
with (root / "metadata/composition/imc_key_audit.tsv").open("w", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(("field", "unique_values", "stage01_exact_key_matches", "rule"))
    writer.writerows(rows)
with (root / "metadata/composition/imc_schema.tsv").open("w", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(("field", "inspection"))
    writer.writerows(schema)
with (root / "metadata/composition/auxiliary_software.tsv").open("w", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(("R_version", "platform", "package", "package_version"))
    writer.writerow(("not_applicable", platform.platform(), "Python", platform.python_version()))
    writer.writerow(("not_applicable", platform.platform(), "Python:h5py", h5py.__version__))

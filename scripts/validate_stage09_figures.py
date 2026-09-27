import sys
import pathlib
import csv
import struct

root = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "external/python"))
import pymupdf

folder = root / "results/stage09_final/pdf_render_audit"
folder.mkdir(exist_ok=True)
rows = []
for file in sorted((root / "figures/final").glob("Figure*.pdf")):
    doc = pymupdf.open(file)
    page = doc[0]
    page.get_pixmap(matrix=pymupdf.Matrix(1.2, 1.2)).save(folder / (file.stem + ".png"))
    outside = sum(1 for x0, y0, x1, y1, *_ in page.get_text("words")
                  if x0 < 0 or y0 < 0 or x1 > page.rect.width or y1 > page.rect.height)
    png = file.with_suffix(".png").read_bytes()
    index = png.index(b"pHYs")
    px, py, unit = struct.unpack(">IIB", png[index+4:index+13])
    dpi = px * .0254
    rows.append(dict(figure=file.stem, pages=len(doc), vector_text_characters=len(page.get_text()),
                     raster_images=len(page.get_images()), out_of_page_words=outside,
                     png_dpi=round(dpi), passed=len(doc) == 1 and outside == 0 and dpi >= 299.9))
with open(root / "results/stage09_final/figure_export_validation.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
assert len(rows) == 6 and all(r["passed"] for r in rows), "Final figure export validation failed"

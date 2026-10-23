"""Retain the saved coverage curves as vectors; clarify denominators without reprocessing raw traces."""
import io,re,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import ContentStream,FloatObject
from reportlab.pdfgen import canvas
from PIL import Image
P=Path(__file__).resolve().parents[2]
ROOT=next(p for p in P.parents if (p/'docs/calibration').is_dir())
SRC=P/'inputs/docs/calibration/figures/sanger_vs_nanopore/figS_sanger_window_original'
OUT=P/'figures/FigureS13_sanger_window'
svg=SRC.with_suffix('.svg').read_text();root=ET.fromstring(svg)
x,y,w,h=map(float,root.attrib['viewBox'].split());W=174/25.4*72;k=W/w;H=h*k+26
svg=re.sub(r'stroke-width:\s*([0-9.]+)',lambda m:f'stroke-width: {max(float(m[1]),.5/k):.6f}',svg)
inner=re.sub(r'^[\s\S]*?<svg[^>]*>','',svg);inner=re.sub(r'</svg>\s*$','',inner)
foot='Coverage: Q ≥ 20 bases per clone. Density: unique events / observed bases within each 25-bp bin.'
OUT.with_suffix('.svg').write_text(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}pt" height="{H}pt" viewBox="0 0 {W} {H}"><rect width="100%" height="100%" fill="white"/><g transform="scale({k})">{inner}</g><text x="10" y="{H-10}" font-family="Arial" font-size="6.5" fill="#4a5562">{foot}</text></svg>')
# Preserve the original vector PDF as well; add the same denominator clarification.
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(W,H));c.setFont('Helvetica',6.5);c.setFillColorRGB(.29,.33,.38)
c.drawString(10,10,foot.replace('≥','>='));c.save();page=PdfReader(buf).pages[0]
reader=PdfReader(SRC.with_suffix('.pdf'));src=reader.pages[0]
kw=W/float(src.mediabox.width)
content=ContentStream(src.get_contents(),reader)
for operands,operator in content.operations:
 if operator==b'w' and operands and 0<float(operands[0])<.5/kw:operands[0]=FloatObject(.5/kw)
src.replace_contents(content)
page.merge_transformed_page(src,Transformation().scale(kw).translate(0,26),over=False)
out=PdfWriter();out.add_page(page);out.write(OUT.with_suffix('.pdf'))
print(OUT)

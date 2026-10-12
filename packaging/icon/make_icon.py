"""Draws the shmtrace app icon (icon.png, 1024 px) and builds icon.icns from it.

Same mark as the web UI's favicon: a navy rounded square with a white N,
plus a read strand carrying three mutation ticks. Run from this folder:
    python make_icon.py
"""
import os, shutil, subprocess
from PIL import Image, ImageDraw

S = 1024
img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
# macOS icon grid: 824 px body centred on a 1024 canvas
m, r = 100, 185
top = Image.new("RGBA", (S, S))
d = ImageDraw.Draw(top)
for y in range(S):  # vertical gradient, #2F5B96 -> #1D3A66
    t = y / S
    c = tuple(int(a + (b - a) * t) for a, b in zip((0x2F, 0x5B, 0x96), (0x1D, 0x3A, 0x66)))
    d.line([(0, y), (S, y)], fill=c + (255,))
mask = Image.new("L", (S, S), 0)
ImageDraw.Draw(mask).rounded_rectangle([m, m, S - m, S - m], r, fill=255)
img.paste(top, (0, 0), mask)

d = ImageDraw.Draw(img)
# the N, from the favicon path M8 23V9h3l6 9V9h3v14h-3l-6-9v9z on a 32 grid
k, ox, oy = 23, 512 - 14 * 23, 250 - 9 * 23
pts = [(8, 23), (8, 9), (11, 9), (17, 18), (17, 9), (20, 9), (20, 23), (17, 23), (11, 14), (11, 23)]
d.polygon([(ox + x * k, oy + y * k) for x, y in pts], fill="white")
# read strand with mutation ticks
y0 = 760
d.rounded_rectangle([250, y0 - 16, 774, y0 + 16], 16, fill=(255, 255, 255, 235))
for x, col in ((360, "#E4572E"), (520, "#F2B134"), (660, "#E4572E")):
    d.ellipse([x - 34, y0 - 34, x + 34, y0 + 34], fill=col, outline="white", width=8)
img.save("icon.png")

iconset = "shmtrace.iconset"
os.makedirs(iconset, exist_ok=True)
for s in (16, 32, 128, 256, 512):
    img.resize((s, s), Image.LANCZOS).save(f"{iconset}/icon_{s}x{s}.png")
    img.resize((2 * s, 2 * s), Image.LANCZOS).save(f"{iconset}/icon_{s}x{s}@2x.png")
subprocess.run(["iconutil", "-c", "icns", iconset, "-o", "icon.icns"], check=True)
shutil.rmtree(iconset)

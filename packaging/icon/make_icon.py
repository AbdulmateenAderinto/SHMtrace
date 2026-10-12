"""Draws the SHMtrace mark and writes every form of it from one definition.

The mark: a sequencing trace (chromatogram peaks) running across a deep-teal tile,
with the tallest peak in coral: a mutation picked out of the trace.

    python make_icon.py        (from this folder)

Writes icon.png (1024 px, macOS icon grid), icon.icns, icon.svg (the same mark as a
vector, used by the README) and favicon.txt (an SVG data URI for the web UI's tab).
"""
import math
import os
import shutil
import subprocess

from PIL import Image, ImageDraw

S = 1024                       # canvas
M, R = 100, 185                # macOS grid: 824-px body, corner radius
TOP, BOT = (0x1F, 0x6F, 0x8B), (0x14, 0x3D, 0x5C)      # teal to deep blue
CORAL, WHITE = (0xFF, 0x6B, 0x57), (255, 255, 255)
# chromatogram peaks along the trace: (x on 0-1, height on 0-1); the 4th is the mutation
PEAKS = [(0.06, 0.42), (0.23, 0.72), (0.40, 0.50), (0.57, 1.00), (0.75, 0.62), (0.93, 0.80)]
MUT = 3
X0, X1, BASE, AMP, SIG = 0.21, 0.79, 0.75, 0.40, 0.055   # trace box on the tile, 0-1 units


def trace(n=240, only=None):
    """Points of the trace (tile units, 0-1); `only` keeps the given peak alone.
    The coral mutation is the trace's own curve around that peak, recoloured."""
    pts = []
    for i in range(n + 1):
        u = i / n
        y = 0.0
        for k, (c, h) in enumerate(PEAKS):
            if only is not None and k != only:
                continue
            y += h * math.exp(-((u - c) ** 2) / (2 * SIG ** 2))
        pts.append((X0 + (X1 - X0) * u, BASE - AMP * min(y, 1.05)))
    return pts


def png(path):
    k = 4                                          # draw large, then downsample: smooth edges
    W = S * k
    tile = Image.new("RGBA", (W, W))
    d = ImageDraw.Draw(tile)
    for y in range(W):
        t = y / W
        c = tuple(int(a + (b - a) * t) for a, b in zip(TOP, BOT, strict=True))
        d.line([(0, y), (W, y)], fill=c + (255,))
    mask = Image.new("L", (W, W), 0)
    ImageDraw.Draw(mask).rounded_rectangle([M * k, M * k, (S - M) * k, (S - M) * k], R * k,
                                           fill=255)
    img = Image.new("RGBA", (W, W), (0, 0, 0, 0))
    img.paste(tile, (0, 0), mask)
    d = ImageDraw.Draw(img)

    def stroke(pts, col, width):
        # a round brush stamped densely along the curve: no seams between segments
        r = width / 2
        for x, y in pts:
            d.ellipse([x * W - r, y * W - r, x * W + r, y * W + r], fill=col)
    lw = int(34 * k)
    stroke(trace(3000), WHITE, lw)
    c, _h = PEAKS[MUT]
    mut = [p for p in trace(3000) if abs((p[0] - X0) / (X1 - X0) - c) < 1.5 * SIG]
    stroke(mut, CORAL, lw)
    img.resize((S, S), Image.LANCZOS).save(path)
    return Image.open(path)


def svg(size=1024, with_grid=True) -> str:
    def path(pts, scale):
        return "M" + " L".join(f"{x * scale:.1f},{y * scale:.1f}" for x, y in pts)
    s = size
    m, r = (M * s / S, R * s / S) if with_grid else (0, 0.22 * s)
    c, _h = PEAKS[MUT]
    mut = [p for p in trace(240) if abs((p[0] - X0) / (X1 - X0) - c) < 1.5 * SIG]
    lw = 34 * s / S
    hexc = "#{:02X}{:02X}{:02X}".format
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {s} {s}">'
            f'<defs><linearGradient id="g" x1="0" y1="0" x2="0" y2="1">'
            f'<stop offset="0" stop-color="{hexc(*TOP)}"/><stop offset="1" stop-color="{hexc(*BOT)}"/>'
            f'</linearGradient></defs>'
            f'<rect x="{m:.1f}" y="{m:.1f}" width="{s - 2 * m:.1f}" height="{s - 2 * m:.1f}" '
            f'rx="{r:.1f}" fill="url(#g)"/>'
            f'<path d="{path(trace(80), s)}" fill="none" stroke="#fff" stroke-width="{lw:.1f}" '
            f'stroke-linecap="round" stroke-linejoin="round"/>'
            f'<path d="{path(mut, s)}" fill="none" stroke="{hexc(*CORAL)}" stroke-width="{lw:.1f}" '
            f'stroke-linecap="round" stroke-linejoin="round"/></svg>')


def main():
    img = png("icon.png")
    with open("icon.svg", "w") as fh:
        fh.write(svg())
    # the tab icon fills its square (no macOS margin), at 64 units
    fav = svg(64, with_grid=False).replace('"', "'").replace("#", "%23").replace("<", "%3C").replace(">", "%3E")
    with open("favicon.txt", "w") as fh:
        fh.write("data:image/svg+xml," + fav)
    iconset = "shmtrace.iconset"
    os.makedirs(iconset, exist_ok=True)
    for s in (16, 32, 128, 256, 512):
        img.resize((s, s), Image.LANCZOS).save(f"{iconset}/icon_{s}x{s}.png")
        img.resize((2 * s, 2 * s), Image.LANCZOS).save(f"{iconset}/icon_{s}x{s}@2x.png")
    subprocess.run(["iconutil", "-c", "icns", iconset, "-o", "icon.icns"], check=True)
    shutil.rmtree(iconset)


if __name__ == "__main__":
    main()

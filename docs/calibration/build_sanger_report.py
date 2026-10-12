"""Build SANGER_VS_NANOPORE.html: the write-up as one page with its figures embedded.

    python docs/calibration/build_sanger_report.py      (needs: pip install markdown)
"""
import base64
import pathlib

import markdown

HERE = pathlib.Path(__file__).resolve().parent
F = HERE / "figures/sanger_vs_nanopore"
FIGS = [  # (section heading the figure closes, file, caption)
    ("## Answer", F / "recommended/fig2_spectrum_concordance.png",
     "Figure 2. Sanger against nanopore, wild-type spectrum (manuscript Figure 3A–E, redrawn from the organised traces). Bars in B and D: pooled proportion with 95% interval; points: each mouse, sized by its mutations."),
    ("## Is the 529-bp window defensible?", F / "figS_sanger_window.png",
     "Figure W. Where the Sanger clones are readable (top) and where the wild-type mutations are (bottom). Blue: the 529-bp window; orange: the two 27-bp extensions of the 556-bp alternatives; dashed: end of the JH4 coding segment."),
    ("## AID dependence on each platform", F / "recommended/fig1_density_wt_vs_ko.png",
     "Figure 1. Wild type against Aicda−/− on each platform (manuscript Figure 2A/B, redrawn from the organised traces)."),
    ("## WT2: the same expanded clone on both platforms", F / "recommended/fig3_agreement_beyond_spectrum.png",
     "Figure 3. Agreement beyond the spectrum: mouse against mouse, position along the intron, hotspot preference, and load per mouse (open: nanopore per clone)."),
    ("## What changed from the manuscript's Sanger numbers", F / "figS_sanger_inclusion.png",
     "Figure S1. What each mouse's Sanger value was built from in the manuscript, and now."),
    (None, F / "recommended/figS_per_mouse_matrices.png",
     "Figure S2. One From→To matrix per mouse and platform (manuscript Supplementary S2/S3)."),
]
CSS = """
:root{--bg:#fbfaf7;--fg:#1d2433;--muted:#5b6475;--rule:#d9dde5;--accent:#2f6ba6;--code:#eef1f6}
@media (prefers-color-scheme:dark){:root{--bg:#14171d;--fg:#e6e9ef;--muted:#a3abba;--rule:#323846;--accent:#7fb0e6;--code:#222733}}
body{background:var(--bg);color:var(--fg);font:16px/1.6 -apple-system,"Helvetica Neue",Arial,sans-serif;margin:0;padding:2.5rem 1rem}
main{max-width:980px;margin:0 auto}
h1{font-size:1.9rem;line-height:1.2;margin:0 0 .5rem}h2{font-size:1.35rem;margin:2.4rem 0 .6rem;padding-top:.6rem;border-top:1px solid var(--rule)}
p,li{max-width:72ch}
table{border-collapse:collapse;margin:1rem 0;font-size:.92rem;display:block;overflow-x:auto}
th,td{border-bottom:1px solid var(--rule);padding:.45rem .7rem;text-align:left;vertical-align:top}th{font-weight:600}
code{background:var(--code);padding:.1rem .3rem;border-radius:4px;font-size:.88em}pre{background:var(--code);padding:.8rem;border-radius:6px;overflow-x:auto}
figure{margin:1.5rem 0;background:#fff;border:1px solid var(--rule);border-radius:8px;padding:.8rem}
figure img{width:100%;height:auto;display:block}figcaption{color:#444;font-size:.88rem;margin-top:.5rem}
a{color:var(--accent)}
"""


def figure(path, cap):
    b = base64.b64encode(path.read_bytes()).decode()
    return f'\n\n<figure><img src="data:image/png;base64,{b}" alt="{cap}"><figcaption>{cap}</figcaption></figure>\n\n'


def main():
    text = (HERE / "SANGER_VS_NANOPORE.md").read_text()
    for heading, path, cap in FIGS:
        if heading is None:
            text += figure(path, cap)
            continue
        i = text.find(heading)
        j = text.find("\n## ", i + len(heading))
        j = len(text) if j < 0 else j
        text = text[:j] + figure(path, cap) + text[j:]
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    html = ('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" '
            'content="width=device-width,initial-scale=1"><title>Sanger vs nanopore, JH4</title>'
            f"<style>{CSS}</style></head><body><main>{body}</main></body></html>")
    (HERE / "SANGER_VS_NANOPORE.html").write_text(html)
    print(HERE / "SANGER_VS_NANOPORE.html")


if __name__ == "__main__":
    main()

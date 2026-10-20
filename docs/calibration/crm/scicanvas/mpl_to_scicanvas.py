"""Turn a matplotlib figure script (crm/figs/*.py) into a SciCanvas builder with one code plot per panel.

    python docs/calibration/crm/scicanvas/mpl_to_scicanvas.py fig3_calibration [...]
    node docs/calibration/crm/scicanvas/build.js docs/calibration/crm/scicanvas/fig3_calibration.js \\
        docs/calibration/crm/figures/Figure3_scicanvas

The figure script runs as usual, but every drawing call it makes on an axes (plot, scatter, bar, text, labels,
limits, ticks, legend, ...) is recorded with its arguments. Each axes becomes one PLOT whose code replays those
calls with the numbers written in, so the panel can be re-run and edited in SciCanvas. The panel letters
(style.letter), the panel titles (set_title) and figure-level text become native SciCanvas text; a figure-level
legend becomes a plot of its own. Each plot's fixed margins put its axes exactly where the matplotlib figure had
them, so panels stay aligned.

Cell Press production rules (text 6-8 pt, strokes 0.5-1.5 pt; build.js checks them) are applied while writing:
any font size under 6 pt or line width under 0.5 pt is raised to the minimum, and every such change is listed.
The builder is regenerated from the analysis outputs, so edits made in SciCanvas belong in the .scifig, and
rerunning this script after new analyses starts from the data again.
"""
import importlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.artist  # noqa: E402
import matplotlib.collections  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.container import Container  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.spines import Spine  # noqa: E402

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "figs"
PAGE_W = 1400                               # SciCanvas page units across the 174-mm figure (build.js)
PRINT_MM = 174
PT_UNIT = PRINT_MM / PAGE_W / 25.4 * 72     # printed points per page unit
MIN_PT, MIN_LW = 6.0, 0.5
SCRIPT_PT = 6.0     # with the style below (6-pt ticks and legends) the worker prints scripts at full size
# The figures' plot style, saved with each converted figure: its smallest sizes are 6 pt, so the plot worker draws
# super- and subscripts at 6 pt (SHRINK_FACTOR = minScriptPt / smallest = 1) instead of shrinking them below it
STYLE = {"fontSize": 6.5, "tickSize": 6, "labelSize": 6.5, "titleSize": 7, "legendSize": 6, "minScriptPt": 6,
         "axisLineWidth": 0.6, "lineWidth": 1.0}
FONT_KW = ("fontsize", "size", "title_fontsize", "labelsize")
LW_KW = ("lw", "linewidth", "linewidths", "markeredgewidth", "mew", "elinewidth", "capthick")
SKIP = {"set_title", "set_position", "get_position", "draw", "apply_aspect"}

_depth = [0]
_on = [False]                  # recording only while the figure script's main() runs
_calls = defaultdict(list)     # id(ax) -> [(name, args, kwargs, ret)]
_fig_calls = []                # (fig, name, args, kwargs, ret)
_ctors = {}                    # id(obj) -> (expr, args, kwargs) for artists the script builds itself
_letters = []                  # (ax, s, dx, dy)
_keep = []


def _top():
    return _on[0] and _depth[0] == 0


def _quiet(cls, name):
    """Calls matplotlib makes inside these (axes creation) are its own, not the script's."""
    orig = getattr(cls, name)

    def w(self, *a, **k):
        _depth[0] += 1
        try:
            return orig(self, *a, **k)
        finally:
            _depth[0] -= 1
    setattr(cls, name, w)


def _wrap(cls, name, sink):
    orig = getattr(cls, name)

    def w(self, *a, **k):
        top = _top()
        _depth[0] += 1
        try:
            r = orig(self, *a, **k)
        finally:
            _depth[0] -= 1
        if top:
            sink(self, name, a, k, r)
        return r
    w.__wrapped__ = orig
    setattr(cls, name, w)


def _ctor(cls, expr):
    orig = cls.__init__

    def init(self, *a, **k):
        top = _top()
        _depth[0] += 1
        try:
            orig(self, *a, **k)
        finally:
            _depth[0] -= 1
        if top:
            _ctors[id(self)] = (expr, a, k)
            _keep.append(self)
    cls.__init__ = init


def install():
    for name in dir(Axes):
        if name.startswith("_") or name.startswith("get_") or name in SKIP | {"set_title"}:
            continue
        if callable(getattr(Axes, name, None)) and name not in ("figure", "axes"):
            _wrap(Axes, name, lambda ax, n, a, k, r: _calls[id(ax)].append((n, a, k, r)))
    _wrap(Axes, "set_title", lambda ax, n, a, k, r: _calls[id(ax)].append((n, a, k, r)))
    for name in ("set_visible", "set_color", "set_linewidth", "set_bounds", "set_edgecolor"):
        _wrap(Spine, name, lambda sp, n, a, k, r: _calls[id(sp.axes)].append(("spine:" + n, (sp,) + a, k, r)))
    for name in ("__init__", "clear", "cla"):
        _quiet(Axes, name)
    for name in ("add_subplot", "add_axes", "subplots", "add_gridspec", "savefig", "draw", "legend", "text"):
        _quiet(Figure, name)
    _ctor(Line2D, "plt.Line2D")
    _ctor(mpatches.Patch, "mpatches.Patch")
    _ctor(mpatches.Rectangle, "mpatches.Rectangle")
    for name in ("legend", "text"):
        _wrap(Figure, name, lambda fig, n, a, k, r: _fig_calls.append((fig, n, a, k, r)))


# ---------------------------------------------------------------- encoding arguments as Python source

class Encoder:
    def __init__(self, ax=None):
        self.ax, self.refs, self.used, self.notes = ax, {}, set(), []

    def num(self, x):
        x = float(x)
        if math.isnan(x):
            return "nan"
        if math.isinf(x):
            return "inf" if x > 0 else "-inf"
        s = f"{x:.6g}"
        return s if any(c in s for c in ".en") else s + ".0" if abs(x) < 1e15 and x != int(x) else s

    def __call__(self, v):
        ax = self.ax
        if id(v) in self.refs:
            self.used.add(self.refs[id(v)].split("[")[0].split(".")[0])
            return self.refs[id(v)]
        if v is None or isinstance(v, (bool, str)):
            return repr(v)
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        if isinstance(v, (float, np.floating)):
            return self.num(v)
        if isinstance(v, np.ma.MaskedArray):
            return self(v.filled(np.nan))
        if isinstance(v, np.ndarray):
            return f"np.array({self(v.tolist())})" if v.ndim else self(v.item())
        if isinstance(v, range):
            return repr(v)
        if isinstance(v, tuple):
            return "(" + ", ".join(self(x) for x in v) + ("," if len(v) == 1 else "") + ")"
        if isinstance(v, list):
            return "[" + ", ".join(self(x) for x in v) + "]"
        if isinstance(v, dict):
            return "{" + ", ".join(f"{self(k)}: {self(x)}" for k, x in v.items()) + "}"
        if ax is not None:
            for expr, t in (("ax.transAxes", ax.transAxes), ("ax.transData", ax.transData),
                            ("ax.get_xaxis_transform()", ax.get_xaxis_transform()),
                            ("ax.get_yaxis_transform()", ax.get_yaxis_transform())):
                if v is t:
                    return expr
        if isinstance(v, mcolors.Colormap):
            cols = [mcolors.to_hex(c, keep_alpha=True) for c in v(np.linspace(0, 1, 64))]
            return f"mcolors.LinearSegmentedColormap.from_list({v.name!r}, {self(cols)})"
        if isinstance(v, mcolors.LogNorm):
            return f"mcolors.LogNorm({self(v.vmin)}, {self(v.vmax)})"
        if isinstance(v, mcolors.Normalize):
            return f"mcolors.Normalize({self(v.vmin)}, {self(v.vmax)})"
        if id(v) in _ctors:
            expr, a, k = _ctors[id(v)]
            return f"{expr}({self.args(a, k)})"
        if isinstance(v, (Line2D, mpatches.Patch, matplotlib.collections.Collection)):
            return self.proxy(v)
        raise TypeError(f"cannot write {type(v).__name__} as code")

    def proxy(self, a):
        """A legend handle drawn elsewhere (another panel), rebuilt from what it looks like."""
        if isinstance(a, Line2D):
            k = dict(color=mcolors.to_hex(a.get_color()), ls=a.get_linestyle(), lw=a.get_linewidth(),
                     marker=a.get_marker(), ms=a.get_markersize(), mfc=a.get_markerfacecolor(),
                     mec=a.get_markeredgecolor())
            return f"plt.Line2D([], [], {', '.join(f'{x}={self(y)}' for x, y in k.items())})"
        if isinstance(a, mpatches.Patch):
            return (f"mpatches.Patch(facecolor={self(mcolors.to_hex(a.get_facecolor(), keep_alpha=True))}, "
                    f"edgecolor={self(mcolors.to_hex(a.get_edgecolor(), keep_alpha=True))})")
        fc = a.get_facecolor()
        ec = a.get_edgecolor()
        return (f"plt.Line2D([], [], ls='', marker='o', ms=4, "
                f"mfc={self(mcolors.to_hex(fc[0]) if len(fc) else 'none')}, "
                f"mec={self(mcolors.to_hex(ec[0]) if len(ec) else 'none')})")

    def args(self, a, k):
        return ", ".join([self(x) for x in a] + [f"{n}={self(x)}" for n, x in k.items()])

    def ref(self, r, name):
        """Name what a call returned, so a later call can pass it on (legend handles)."""
        self.refs[id(r)] = name
        if isinstance(r, (list, tuple, Container)):
            for i, x in enumerate(r):
                self.refs.setdefault(id(x), f"{name}[{i}]")


def _scripts(a, k):
    """Whether a call's text carries mathtext super- or subscripts."""
    texts = [x for x in list(a) + list(k.values()) if isinstance(x, str)]
    texts += [y for x in a if isinstance(x, (list, tuple)) for y in x if isinstance(y, str)]
    return any("$" in t and ("^" in t or "_" in t) for t in texts)


def comply(name, kwargs, notes, where, args=()):
    """Cell Press minimums: fonts >= 6 pt (scripts too), lines >= 0.5 pt."""
    k = dict(kwargs)
    floor = SCRIPT_PT if _scripts(args, k) else MIN_PT
    for key in FONT_KW:
        if isinstance(k.get(key), (int, float)) and k[key] < floor:
            notes.append(f"{where}: {name} font {k[key]} -> {floor} pt")
            k[key] = floor
    for key in LW_KW:
        v = k.get(key)
        if isinstance(v, (int, float)) and 0 < v < MIN_LW:
            notes.append(f"{where}: {name} {key} {v} -> {MIN_LW} pt")
            k[key] = MIN_LW
    return k


# ---------------------------------------------------------------- layout

def _page(fig):
    w_in, h_in = fig.get_size_inches()
    s = PAGE_W / (w_in * 25.4)                     # page units per mm
    return s, h_in * 25.4 * s, fig.dpi


def _box(bbox, fig):
    """A display-space bbox as page (x, y, w, h), y down."""
    s, H, dpi = _page(fig)
    mm = 25.4 / dpi
    return (bbox.x0 * mm * s, H - bbox.y1 * mm * s, bbox.width * mm * s, bbox.height * mm * s)


def markup(s):
    """matplotlib mathtext in a title -> SciCanvas rich text (J$_H$4 -> J_{H}4)."""
    def sub(m):
        t = m.group(1)
        t = re.sub(r"_\{([^}]*)\}|_(\w)", lambda q: "_{" + (q.group(1) or q.group(2)) + "}", t)
        t = re.sub(r"\^\{([^}]*)\}|\^(\w)", lambda q: "^{" + (q.group(1) or q.group(2)) + "}", t)
        t = re.sub(r"\\mathit\{([^}]*)\}", r"{i|\1}", t)
        return t.replace("\\mu", "μ").replace("\\", "")
    return re.sub(r"\$([^$]*)\$", sub, s)


def js(s):
    return json.dumps(s, ensure_ascii=False)


# ---------------------------------------------------------------- one figure

def convert(stem, figure_name, out_js):
    _calls.clear(), _fig_calls.clear(), _ctors.clear(), _letters.clear()
    sys.path.insert(0, str(FIGS))
    import style
    rc = {}
    real_update = plt.rcParams.update
    style_setup = style.setup

    def capture_setup():
        def grab(d=None, **k):
            rc.update(d or {}, **k)
            real_update(d or {}, **k)
        plt.rcParams.update = grab
        try:
            style_setup()
        finally:
            plt.rcParams.update = real_update
    mod = importlib.import_module(stem)
    captured = {}
    mod.setup = capture_setup
    mod.save = lambda fig, name: captured.setdefault("fig", fig)
    mod.letter = lambda ax, s, dx=-0.18, dy=1.06: _letters.append((ax, s, dx, dy))
    _on[0] = True
    try:
        mod.main()
    finally:
        _on[0] = False
    fig = captured["fig"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    s, H, dpi = _page(fig)
    notes = []
    rc_ok = {k: v for k, v in rc.items() if k not in ("pdf.fonttype", "ps.fonttype", "svg.fonttype")}
    # what the original drew with matplotlib's defaults, where the plot worker's style would differ
    for k in ("xtick.minor.width", "ytick.minor.width", "patch.linewidth"):
        rc_ok.setdefault(k, max(MIN_LW, matplotlib.rcParamsDefault[k]))
    for k in list(rc_ok):
        if k.endswith(("labelsize", "fontsize", "titlesize")) or k == "font.size":
            if isinstance(rc_ok[k], (int, float)) and rc_ok[k] < MIN_PT:
                notes.append(f"style: {k} {rc_ok[k]} -> {MIN_PT} pt")
                rc_ok[k] = MIN_PT
    out = [f"() => {{\n  // {figure_name}, rebuilt from figs/{stem}.py by mpl_to_scicanvas.py: one code plot per panel,\n"
           f"  // the panel letters and titles as SciCanvas text. Regenerate with that script; build with build.js.\n"
           f"  const W = {PAGE_W}, H = {math.ceil(H)}, o = [];\n"
           f"  STYLE({json.dumps(STYLE)});\n"
           "  const F = 'arial', INK = '#1d2733';\n"
           "  // text kept on the page: a title that ran past the edge in matplotlib is moved in, not cropped\n"
           "  const T = (text, x, y, e = {}) => { const t = Make.text(text, x, y, { color: INK, family: F, ...e }); "
           "t.x = Math.max(0, Math.min(t.x, W - t.w)); o.push(t); return t; };\n"
           "  const B = (text, x, bottom, e = {}) => { const t = T(text, x, 0, e); t.y = bottom - t.h; return t; };\n"]
    header = ("import matplotlib.colors as mcolors\nimport matplotlib.patches as mpatches\n"
              "nan, inf = float('nan'), float('inf')\n"
              f"plt.rcParams.update({rc_ok!r})   # the figures' style (figs/style.py)\n")
    letters = {id(ax): (t, dx, dy) for ax, t, dx, dy in _letters}
    for i, ax in enumerate(fig.axes):
        calls = _calls.get(id(ax), [])
        title = next(((a, k) for n, a, k, _r in calls if n == "set_title"), None)
        letter = letters.get(id(ax))
        name = (letter[0] if letter else f"panel {i + 1}") + (f": {markup(title[0][0])}" if title else "")
        enc = Encoder(ax)
        lines = []
        grow = 1.0                  # how much this panel's legend text grew to reach 6 pt
        for n, a, k, r in calls:
            if n == "set_title" or n.startswith("set_position"):
                continue
            seen = a
            if n == "legend":       # a legend's labels come from earlier calls' label= arguments
                seen = tuple(a) + (ax.get_legend_handles_labels()[1],)
                if "fontsize" not in k and _scripts(seen, {}):
                    k = {**k, "fontsize": plt.rcParams["legend.fontsize"]}
            before = k.get("fontsize")
            k = comply(n, k, notes, name, seen)
            if n == "legend" and k.get("fontsize") != before:
                old = before if isinstance(before, (int, float)) else plt.rcParams["legend.fontsize"]
                grow = max(grow, k["fontsize"] / old)
                # a legend raised to 6 pt grows wider than the space it had: shorter handles win it back
                k.setdefault("handlelength", 1.5)
                k["handletextpad"] = min(k.get("handletextpad", 0.5), 0.5)
            if n in ("set_xlabel", "set_ylabel") and _scripts(a, k) and "fontsize" not in k \
                    and plt.rcParams["axes.labelsize"] < SCRIPT_PT:
                k["fontsize"] = SCRIPT_PT
                notes.append(f"{name}: {n} with a superscript -> {SCRIPT_PT} pt")
            if n in ("set_xscale", "set_yscale") and a and a[0] == "log" and SCRIPT_PT > MIN_PT:
                # log tick labels are 10^n: their exponents need the larger size
                lines.append(f"ax.{n}({enc.args(a, k)})")
                lines.append(f"ax.tick_params(axis={n[4]!r}, labelsize={SCRIPT_PT})")
                notes.append(f"{name}: log {n[4]} tick labels -> {SCRIPT_PT} pt (their exponents print at 6)")
                continue
            if n.startswith("spine:"):
                sp, rest = a[0], a[1:]
                key = next(x for x, y in ax.spines.items() if y is sp)
                lines.append(f"ax.spines[{key!r}].{n[6:]}({enc.args(rest, k)})")
                continue
            try:
                call = f"ax.{n}({enc.args(a, k)})"
            except TypeError as e:
                raise TypeError(f"{name}: ax.{n}: {e}") from e
            var = f"r{len(lines)}"
            if isinstance(r, (list, Container, matplotlib.artist.Artist)):
                enc.ref(r, var)
            lines.append((var, call))
        code_lines = []
        for x in lines:
            if isinstance(x, tuple):
                v, c = x
                code_lines.append(f"{v} = {c}" if v in enc.used else c)
            else:
                code_lines.append(x)
        # the panel without its title (drawn as SciCanvas text), as it sits on the page
        titles = (ax.title, ax._left_title, ax._right_title)
        shown = [t.get_visible() for t in titles]
        for t in titles:
            t.set_visible(False)
        tight = ax.get_tightbbox(renderer)
        for t, v in zip(titles, shown, strict=True):
            t.set_visible(v)
        px, py, pw, ph = _box(tight, fig)
        # room around the panel's own text: 3 mm across (text raised to 6 pt grows sideways), 1.5 mm up and down
        hpad, vpad = 3.0 * s, 1.5 * s
        x0 = max(0, px - hpad)
        y0 = max(0, py - vpad)
        px, py, pw, ph = x0, y0, min(PAGE_W - 6 - x0, pw + (px - x0) + hpad), min(H - 6 - y0, ph + (py - y0) + vpad)
        axx, axy, axw, axh = _box(ax.get_window_extent(renderer), fig)
        leg = ax.get_legend()
        if leg is not None and grow > 1:
            lx, ly, lw_, lh = _box(leg.get_window_extent(renderer), fig)
            extra = lw_ * (grow - 1) + 1.0 * s
            if lx + lw_ > axx + axw:          # a legend beside the axes: give its larger text room
                pw = min(PAGE_W - px, pw + extra)
        margins = {"left": (axx - px) / s, "right": (px + pw - axx - axw) / s,
                   "top": (axy - py) / s, "bottom": (py + ph - axy - axh) / s}
        code = header + "fig, ax = plt.subplots()\n" + "\n".join(code_lines) + "\n"
        out.append(f"  // {name}\n  o.push(PLOT({{ name: {js(name)}, x: {px:.1f}, y: {py:.1f}, w: {pw:.1f}, h: {ph:.1f}, "
                   f"margins: {json.dumps({k: round(v, 2) for k, v in margins.items()})}, code: {js(code)} }}));\n")
        if title:
            # where matplotlib puts a left title: at the axes' left edge, `pad` points above its top
            (text,), tk = title[0], title[1]
            size = max(MIN_PT, tk.get("fontsize", plt.rcParams["axes.titlesize"]))
            pad_pt = tk.get("pad", plt.rcParams["axes.titlepad"])
            bottom = axy - pad_pt / 72 * 25.4 * s
            out.append(f"  B({js(markup(text))}, {axx:.1f}, {bottom:.1f}, "
                       f"{{ fontSize: {size / PT_UNIT:.2f}, bold: true }});\n")
        if letter:
            lt, dx, dy = letter
            disp = ax.transAxes.transform((dx, dy))
            bx, by = disp[0] * 25.4 / dpi * s, H - disp[1] * 25.4 / dpi * s
            out.append(f"  B({js(lt)}, {bx:.1f}, {by:.1f}, {{ fontSize: {8 / PT_UNIT:.2f}, bold: true }});\n")
    for fig_, n, a, k, r in _fig_calls:
        if fig_ is not fig:
            continue
        if n == "legend":
            enc = Encoder(None)
            labels = [h.get_label() for h in (k.get("handles") or [])] + list(k.get("labels") or [])
            if "fontsize" not in k and _scripts((labels,), {}):
                k = {**k, "fontsize": plt.rcParams["legend.fontsize"]}
            k = comply("legend", k, notes, "figure legend", (labels,))
            k = {x: y for x, y in k.items() if x not in ("loc", "bbox_to_anchor")}
            box = list(_box(r.get_window_extent(renderer), fig))
            slack = 4 * s                       # a little room each side: the legend is centred in its plot
            box[0], box[2] = max(0, box[0] - slack), min(PAGE_W - max(0, box[0] - slack), box[2] + 2 * slack)
            code = (header + "fig, ax = plt.subplots()\nax.set_axis_off()\n"
                    f"ax.legend({enc.args(a, k)}, loc='center')\n")
            out.append(f"  o.push(PLOT({{ name: 'legend', x: {box[0]:.1f}, y: {box[1]:.1f}, w: {box[2]:.1f}, "
                       f"h: {box[3]:.1f}, margins: {{ left: 0, right: 0, top: 0, bottom: 0 }}, code: {js(code)} }}));\n")
        elif n == "text":
            tb = _box(r.get_window_extent(renderer), fig)
            out.append(f"  T({js(markup(r.get_text()))}, {tb[0]:.1f}, {tb[1]:.1f}, "
                       f"{{ fontSize: {max(MIN_PT, r.get_fontsize()) / PT_UNIT:.2f} }});\n")
    out.append(f"  return {{ id: {js('crm-' + stem)}, name: {js(figure_name)}, width: W, height: H, "
               f"background: '#ffffff', objects: o, alt: {js(figure_name + ', rebuilt from its matplotlib script.')} }};\n}}\n")
    Path(out_js).write_text("".join(out))
    plt.close(fig)
    sys.modules.pop(stem, None)
    return notes


FIGURES = {"fig3_calibration": "Figure 3", "fig4_controls": "Figure 4", "fig7_platforms": "Figure 7",
           "figS3_estimates": "Figure S3", "figS5_depth": "Figure S5"}

if __name__ == "__main__":
    install()
    for stem in sys.argv[1:] or FIGURES:
        notes = convert(stem, FIGURES.get(stem, stem), HERE / f"{stem}.js")
        print(f"{stem}.js written; {len(notes)} Cell Press adjustment(s)")
        for n in notes:
            print("   ", n)

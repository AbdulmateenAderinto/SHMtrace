"""Renumber the supplementary figures and tables into the order the article first cites them (Cell Press rule).

    python crm/manuscript/renumber_supp.py check          # print the citation order and the mapping it implies
    python crm/manuscript/renumber_supp.py apply          # rewrite every citation in the sources at once

The order is read from the built article (crm/build/SHMtrace_CellReportsMethods.docx, through macOS textutil), up
to the main figure legends, so the Summary, Results, Discussion and STAR Methods all count. `apply` rewrites
"Figure(s) S…" and "Table(s) S…" citations, legend titles included, in one pass, so S19 -> S4 and S4 -> S20 cannot
collide; a range whose members stop being consecutive is written out as a list. Citations split across Python
string literals ("Figures S8 " "and S9") are followed through the split. File names, script names and placeholder
keys are not renumbered: build_crm.py names each legend's image explicitly. Rebuild and run `check` again after
`apply`: the mapping should then be the identity.
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CRM = HERE.parent
ARTICLE = CRM / "build" / "SHMtrace_CellReportsMethods.docx"
SOURCES = [HERE / f for f in ("crm_text.py", "crm_methods.py", "build_crm.py", "crm_numbers.py", "lab_details.py",
                              "assemble.py")] + [CRM / "make_tables.py", CRM / "scicanvas/figS20_platforms.js",
                                                 CRM / "scicanvas/figS21_overview.js"]
DASH = r"(?:–|\\u2013|-)"
JOIN = r'(?:"\s*\n\s*f?")?'                         # a citation running on into the next string literal
ITEM = rf"S?\d+[A-Z]*(?:{JOIN}{DASH}{JOIN}S?\d+[A-Z]*)?"
SEP = rf"{JOIN}(?:,\s+and\s+|,\s+|\s+and\s+|;\s+){JOIN}"
CITE = re.compile(rf"\b(Figures?|Tables?)(\s+){JOIN}({ITEM}(?:{SEP}{ITEM})*)")


def items(text):
    """[(supplementary?, number, letters, range_end or None)] of a citation's item list."""
    out = []
    for m in re.finditer(rf"(S?)(\d+)([A-Z]*)(?:{JOIN}{DASH}{JOIN}(S?)(\d+)([A-Z]*))?", text):
        out.append((m.group(1) == "S", int(m.group(2)), m.group(3), int(m.group(5)) if m.group(5) else None))
    return out


def order():
    text = subprocess.run(["textutil", "-convert", "txt", "-stdout", str(ARTICLE)], capture_output=True,
                          text=True, check=True).stdout
    text = re.sub(r"\s+", " ", text)
    end = text.find("Figure 1. SHMtrace")
    body = text[:end] if end > 0 else text
    seen = {"Figure": [], "Table": []}
    for m in CITE.finditer(body):
        kind = "Figure" if m.group(1).startswith("Fig") else "Table"
        for sup, n, _l, hi in items(m.group(3)):
            if sup:
                for k in range(n, (hi or n) + 1):
                    if k not in seen[kind]:
                        seen[kind].append(k)
    return seen


def mapping(seen):
    return {kind: {old: new for new, old in enumerate(seq, 1)} for kind, seq in seen.items()}


def render(kind_word, entries):
    """Write mapped (number, letters) entries as a citation list, consecutive whole figures as ranges."""
    entries = sorted(set(entries))
    parts, i = [], 0
    while i < len(entries):
        n, letters = entries[i]
        j = i
        while (not letters and j + 1 < len(entries) and entries[j + 1][0] == entries[j][0] + 1
               and not entries[j + 1][1]):
            j += 1
        parts.append(f"S{n}{letters}" if j - i < 2 else f"S{n}–S{entries[j][0]}")
        if 0 < j - i < 2:                      # two consecutive figures read better as "S8 and S9"
            parts.append(f"S{entries[j][0]}")
        i = j + 1
    plural = kind_word + ("s" if len(parts) > 1 or "–" in "".join(parts) else "")
    body = parts[0] if len(parts) == 1 else ", ".join(parts[:-1]) + (" and " if len(parts) == 2 else ", and ") + parts[-1]
    return plural, body


def rewrite(src, maps):
    def one(m):
        word = "Figure" if m.group(1).startswith("Fig") else "Table"
        its = items(m.group(3))
        if not any(sup for sup, *_ in its):
            return m.group(0)
        if any(not sup for sup, *_ in its):    # a list mixing main and supplementary items: renumber in place
            def num(mm):
                return f"S{maps[word].get(int(mm.group(1)), int(mm.group(1)))}"
            return re.sub(r"S(\d+)", num, m.group(0))
        mp, entries = maps[word], []
        for _sup, n, letters, hi in its:
            if hi is None:
                entries.append((mp.get(n, n), letters))
            else:
                entries += [(mp.get(k, k), "") for k in range(n, hi + 1)]
        unchanged = all(mp.get(n, n) == n and (hi is None or all(mp.get(k, k) == k for k in range(n, hi + 1)))
                        for _s, n, _l, hi in its)
        if unchanged:
            return m.group(0)
        plural, body = render(word, entries)
        # keep the original noun's number when it was already plural (lists), else follow the new list
        noun = m.group(1) if m.group(1).endswith("s") else plural
        return f"{noun}{m.group(2)}{body}"
    return CITE.sub(one, src)


def main(cmd):
    seen = order()
    maps = mapping(seen)
    for kind in ("Figure", "Table"):
        moved = {o: n for o, n in maps[kind].items() if o != n}
        print(f"{kind}s S, first-citation order: {seen[kind]}")
        print(f"  renumbering: {moved or 'none needed'}")
    if cmd != "apply":
        return
    for f in SOURCES:
        if not f.exists():
            continue
        s = f.read_text()
        new = rewrite(s, maps)
        if new != s:
            f.write_text(new)
            print("rewrote", f.relative_to(CRM))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "check")

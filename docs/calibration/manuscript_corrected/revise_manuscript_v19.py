"""Ver19: Ver18 plus mutation trees (Supplementary Figure S10).

    python revise_manuscript_v19.py VER18.docx VER19.docx      (needs python-docx)

Adds a Methods paragraph (how trees and the four-gamete test are made), a Results paragraph
in Section 3.2, and Supplementary Figure S10 (docs/calibration/mutation_tree_figure.py,
mutation_trees_wt3.png) after S9. New text is highlighted yellow. SOURCE is never
modified, and the script stops if an anchor is not found.
"""
import copy
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_COLOR_INDEX
from docx.shared import Inches
from docx.text.paragraph import Paragraph

from revise_manuscript_v18 import restyle

HERE = Path(__file__).resolve().parent
TREES = HERE.parent / "figures/mutation_trees/mutation_trees_wt3.png"

METHODS = (
    "Mutation trees. In the full-length analysis, reads were grouped into clones by their V(D)J junction "
    "(the 20 bases 5′ of the JH4 coding segment, at most 3 edits apart). Within a clone, every call carried by at "
    "least two reads was listed with the reads carrying it. Calls carried by nearly the same reads (Jaccard index "
    "≥ 0.9) formed one node, because the reads cannot order them, and each node's parent was the smallest node whose "
    "reads included at least 90% of its own: the mutation that arose before it on its lineage. A read fitted the "
    "tree when every node it carried lay on one path from the root. Tree-likeness was also measured with the "
    "four-gamete test: two calls, each carried by at least three reads, are incompatible with any single tree when "
    "reads carrying both, either one alone and neither are all present (each in at least two reads).")

RESULTS = (
    "Because each read spans the amplicon, the reads of one B-cell clone can also be ordered into a mutation tree "
    "(Supplementary Figure S10). The largest Wt3 clone (84 reads sharing a V(D)J junction) carried four founder "
    "mutations in every mutated read and then divided into five sub-lineages, up to 16 mutations deep; 87% of its "
    "mutated reads fitted a single tree, and 0.9% of its mutation pairs failed the four-gamete test. Every other "
    "Wt clone with at least 10 reads had no incompatible pairs. In the other Wt mice, the dominant Wt2 clone was unmutated, and each Wt4 clone carried one "
    "block of shared mutations without sub-branches, also at 911 reads. Without UMIs, recombination between "
    "molecules during PCR can create incompatible pairs, so these trees are read-level reconstructions rather than "
    "confirmed genealogies.")

LEGEND = (
    "Mutation tree of a Wt3 JH4 clone from full-length reads. (A) The largest Wt3 clone (84 reads sharing a V(D)J "
    "junction; full-length analysis, positions 90–1,230, 500 reads). Square, the clone's unmutated ancestor; ticks, "
    "mutations, coloured G:C transition, G:C transversion or A:T; circles, nodes (mutations carried by the same "
    "reads), whose area and number give the reads carrying them; horizontal distance, mutations accumulated since "
    "the ancestor. Mutations on one branch cannot be ordered by the reads and are drawn in position order. 87% of "
    "the clone's mutated reads fit the tree, and 10 of 1,128 mutation pairs fail the four-gamete test. (B) The same "
    "clone's reads (rows) by mutations (columns, in tree order): each sub-lineage forms a block that keeps every "
    "earlier mutation. Without UMIs, PCR recombination between molecules can produce reads that fit no branch, so "
    "the tree is a read-level reconstruction, not a confirmed genealogy.")


def insert_after(anchor, text):
    """A new paragraph after `anchor`, in its style, highlighted yellow."""
    p = copy.deepcopy(anchor._p)
    anchor._p.addnext(p)
    np_ = Paragraph(p, anchor._parent)
    for r in list(np_.runs)[1:]:
        p.remove(r._r)
    for child in list(p):            # drop hyperlinks and bookmarks copied from the anchor
        if child.tag.split("}")[1] not in ("pPr", "r"):
            p.remove(child)
    run = np_.runs[0]
    run.text = text
    run.italic = run.bold = None
    run.font.subscript = run.font.superscript = None
    run.font.highlight_color = WD_COLOR_INDEX.YELLOW
    return np_


def main(src, out):
    doc = Document(src)

    def find(start):
        hits = [p for p in doc.paragraphs if p.text.startswith(start)]
        if len(hits) != 1:
            raise RuntimeError(f"{len(hits)} paragraphs start with {start!r}")
        return hits[0]

    restyle(insert_after(find("For the separate full-length ONT analysis"), METHODS))
    restyle(insert_after(find("Because ONT reads span the complete amplicon"), RESULTS))

    # Supplementary Figure S10 after S9, in S9's style: picture, then legend
    s9 = find("Supplementary Figure S9.")
    pic = copy.deepcopy(s9._p)
    for child in list(pic):
        if child.tag.split("}")[1] != "pPr":
            pic.remove(child)
    s9._p.addnext(pic)
    Paragraph(pic, s9._parent).add_run().add_picture(str(TREES), width=Inches(6.4))
    leg = copy.deepcopy(s9._p)
    pic.addnext(leg)
    lp = Paragraph(leg, s9._parent)
    bold, body = lp.runs[0], lp.runs[1]
    for r in lp.runs[2:]:
        leg.remove(r._r)
    for child in list(leg):
        if child.tag.split("}")[1] not in ("pPr", "r"):
            leg.remove(child)
    bold.text, body.text = "Supplementary Figure S10. ", LEGEND
    for r in (bold, body):
        r.font.highlight_color = WD_COLOR_INDEX.YELLOW
    restyle(lp)
    doc.core_properties.comments = "Ver19: Ver18 plus mutation trees (Supplementary Figure S10)"
    doc.save(out)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

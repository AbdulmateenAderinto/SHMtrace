"""Mutation trees: the order in which a clone's mutations arose, read from long reads.

Every read of a B-cell clone descends from the clone's founder, so a mutation that arose
early is carried by every read of the branch below it, and a later mutation by a subset
of those reads. Long reads keep all of a molecule's mutations together, so the nesting
can be read directly:

1. Within one clone (reads sharing a V(D)J junction), each called mutation is listed
   with the reads that carry it. Mutations seen in fewer than `min_reads` reads are left
   out of the tree: one read cannot place a mutation, and it may be an error.
2. Mutations carried by (nearly) the same reads arose on the same branch and cannot be
   ordered from these reads; they form one node (Jaccard >= `merge`).
3. A node's parent is the smallest node whose reads contain at least `contain` of its
   own reads: the mutation that arose before it on its lineage. A node contained by no
   other hangs from the root, the clone's unmutated ancestor.

The tolerance (`contain` < 1) absorbs sequencing errors and dropouts at a site. PCR
recombination between molecules would produce reads that fit no branch; such reads
show up as a low `fit` (share of carrier reads consistent with the tree).

    tree = build({read: {(pos, alt), ...}, ...})
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field


@dataclass
class Node:
    id: int
    mutations: list = field(default_factory=list)   # [(pos, alt), ...], 1-based positions
    reads: frozenset = frozenset()                  # reads carrying the node's mutations
    parent: int | None = None                       # None: the root (unmutated ancestor)
    children: list = field(default_factory=list)


def build(read_mutations: dict, min_reads: int = 2, merge: float = 0.9,
          contain: float = 0.9) -> dict:
    """Mutation tree of one clone. `read_mutations`: read name -> set of (pos, alt)."""
    carriers = defaultdict(set)
    for r, muts in read_mutations.items():
        for m in muts:
            carriers[m].add(r)
    kept = sorted((m for m, c in carriers.items() if len(c) >= min_reads),
                  key=lambda m: (-len(carriers[m]), m))
    nodes: list[Node] = []
    for m in kept:
        c = carriers[m]
        for n in nodes:
            inter = len(c & n.reads)
            if inter / len(c | n.reads) >= merge:
                n.mutations.append(m)
                break
        else:
            nodes.append(Node(len(nodes) + 1, [m], frozenset(c)))
    for i, n in enumerate(nodes):
        best = None
        for a in nodes[:i]:
            if len(a.reads) >= len(n.reads) and len(n.reads & a.reads) >= contain * len(n.reads):
                if best is None or len(a.reads) < len(best.reads):
                    best = a
        if best is not None:
            n.parent = best.id
            best.children.append(n.id)
    by_id = {n.id: n for n in nodes}
    # a read sits at the deepest node it carries; reads with no kept mutation sit at the root
    at = Counter()
    fit = []
    for r, muts in read_mutations.items():
        mine = [n for n in nodes if r in n.reads]
        if not mine:
            at[0] += 1
            continue
        deepest = max(mine, key=lambda n: _depth(n, by_id))
        at[deepest.id] += 1
        # consistent: every node the read carries is on the deepest node's path to the root
        path = set(_path(deepest, by_id))
        fit.append(all(n.id in path for n in mine))
    for n in nodes:
        n.mutations.sort()
    return {"nodes": nodes, "reads": len(read_mutations), "reads_at": dict(at),
            "fit": round(sum(fit) / len(fit), 3) if fit else None,
            "left_out": sum(1 for c in carriers.values() if len(c) < min_reads)}


def _path(n: Node, by_id: dict) -> list[int]:
    out = [n.id]
    while n.parent is not None:
        n = by_id[n.parent]
        out.append(n.id)
    return out


def _depth(n: Node, by_id: dict) -> int:
    """Mutations between the root and the end of this node's branch."""
    return sum(len(by_id[i].mutations) for i in _path(n, by_id))


def depth(tree: dict, node_id: int) -> int:
    by_id = {n.id: n for n in tree["nodes"]}
    return _depth(by_id[node_id], by_id)


def clone_reads(aligned, read_rows: list, variants: list, p, max_edits: int | None = None) -> dict:
    """Group a sample's reads by V(D)J junction and give each read its called mutations.

    `read_rows` is pipeline.per_read_table's output for the same reads (mutation
    positions); `variants` the sample's calls (the allele at each position).
    Returns {clone id: {read: {(pos, alt), ...}}}, largest clone first.
    """
    from . import clones
    alt_at = defaultdict(set)
    for v in variants:
        alt_at[v.pos].add(v.alt)
    j0 = p.junction_ref_start - 1
    keys = {}
    for r in aligned:
        k = clones.junction(r, j0, p.junction_length)
        if k:
            keys[r.qname] = k
    clone_of = clones.cluster(keys, p.junction_max_edits if max_edits is None else max_edits)
    base_at = {r.qname: {rpos + 1: b for rpos, b, _q in r.pairs} for r in aligned if r.qname in clone_of}
    out = defaultdict(dict)
    for row in read_rows:
        name = row["read"]
        if name not in clone_of:
            continue
        muts = set()
        for pos in row["positions"]:
            b = base_at[name].get(pos)
            if b in alt_at.get(pos, ()):
                muts.add((pos, b))
        out[clone_of[name]][name] = muts
    return dict(sorted(out.items(), key=lambda kv: -len(kv[1])))


__all__ = ["Node", "build", "depth", "clone_reads"]

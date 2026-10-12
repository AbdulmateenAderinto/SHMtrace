"""Read input for every supported platform, plus Phred/expected-error math.

`read_fastq(path)` is the one entry point; it yields (name, seq, qual) tuples
whatever the input is:

  * FASTQ (.fastq/.fq), any platform
  * FASTA (.fasta/.fa/.fna/.fas), which has no qualities: every base is given
    the constant Phred FASTA_DEFAULT_Q (30, i.e. 1 error in 1,000)
  * Sanger traces (.ab1/.abi): the basecalls and per-base Phred values stored
    in the ABIF file (a built-in parser, no Biopython needed)
  * plain-text Sanger sequence (.seq): one read per file, named after the file,
    with FASTA_DEFAULT_Q qualities (a .seq that is really FASTA is read as FASTA)
  * GenBank (.gb/.gbk), as sequencing providers and Benchling export Sanger
    results: each record's sequence, with FASTA_DEFAULT_Q qualities. IUPAC
    ambiguity codes are kept, so they are never counted as mutations
  * BAM/SAM/CRAM, when pysam is installed (`pip install nanoshm[bam]`)
  * gzip of any of the text formats (add .gz)
  * Illumina paired-end reads as "R1.fastq.gz+R2.fastq.gz": the two files are
    merged pair by pair into single reads (merge_pairs)
  * a directory: every readable file inside it is one sample (for example a
    folder of Sanger clones). Files whose names differ only in _R1/_R2 (also
    _1/_2 before the extension) are merged as Illumina pairs.

Files of unknown extension are sniffed by content (ABIF magic, '@' FASTQ,
'>' FASTA). In a directory, files that are none of these are skipped.
"""
from __future__ import annotations

import gzip
import os
import re
import struct
from collections.abc import Iterator

# Phred assigned to every base of a read that comes without qualities
# (FASTA, .seq). Q30 = one error per 1,000 bases: passes the default EE gate
# and the Sanger Q20 calling threshold, so such reads are trusted as given.
FASTA_DEFAULT_Q = 30

FASTQ_EXT = (".fastq", ".fq")
FASTA_EXT = (".fasta", ".fa", ".fna", ".fas", ".fsa")
ABIF_EXT = (".ab1", ".abi", ".ab!")
SEQ_EXT = (".seq",)
GENBANK_EXT = (".gb", ".gbk", ".genbank")
ALIGN_EXT = (".bam", ".sam", ".cram")
READ_EXTENSIONS = FASTQ_EXT + FASTA_EXT + ABIF_EXT + SEQ_EXT + GENBANK_EXT + ALIGN_EXT

PAIR_SEPARATOR = "+"


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path: str) -> list[tuple[str, str]]:
    recs, name, seq = [], None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = (line[1:].split() or ["unnamed"])[0], []
            elif line:
                seq.append(line.strip().upper())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def is_alignment_file(path: str) -> bool:
    return str(path).lower().endswith(ALIGN_EXT)


def read_bam(path: str) -> Iterator[tuple[str, str, str]]:
    """Reads from a BAM/SAM as (name, seq, qual) in their sequenced orientation.

    Aligners store reverse-strand reads reverse-complemented; this undoes that,
    so nanoshm realigns the reads exactly as they came off the sequencer. A read
    with several primary records (Geneious writes split alignments this way)
    is returned once, from its longest record. Hard-clipped records are
    skipped because the clipped bases are missing.
    """
    try:
        import pysam
    except ImportError as e:                       # pragma: no cover
        raise ImportError("reading BAM/SAM needs pysam: pip install 'nanoshm[bam]'") from e
    best: dict[str, tuple[str, str]] = {}
    order: list[str] = []
    with pysam.AlignmentFile(path, check_sq=False) as fh:
        for a in fh.fetch(until_eof=True):
            if a.is_secondary or a.is_supplementary or not a.query_sequence:
                continue
            if a.cigartuples and any(op == 5 for op, _ in a.cigartuples):
                continue
            seq = a.query_sequence.upper()
            q = a.query_qualities
            qual = "".join(chr(x + 33) for x in q) if q is not None else "I" * len(seq)
            if a.is_reverse:
                seq, qual = seq.translate(_COMP)[::-1], qual[::-1]
            name = a.query_name
            if name not in best:
                order.append(name)
            if name not in best or len(seq) > len(best[name][0]):
                best[name] = (seq, qual)
    for name in order:
        yield name, best[name][0], best[name][1]


# ---------------------------------------------------------------- formats --

def _strip_gz(path: str) -> str:
    p = path.lower()
    return p[:-3] if p.endswith(".gz") else p


def _stem(path: str) -> str:
    base = os.path.basename(path.rstrip("/"))
    if base.lower().endswith(".gz"):
        base = base[:-3]
    return os.path.splitext(base)[0] or base


def sniff_format(path: str) -> str | None:
    """'fastq', 'fasta', 'abif', 'seq', 'genbank', 'bam' or None (not a read file)."""
    low = _strip_gz(path)
    if low.endswith(ALIGN_EXT):
        return "bam"
    if low.endswith(ABIF_EXT):
        return "abif"
    if low.endswith(FASTQ_EXT):
        return "fastq"
    if low.endswith(FASTA_EXT):
        return "fasta"
    if low.endswith(SEQ_EXT):
        return "seq"
    if low.endswith(GENBANK_EXT):
        return "genbank"
    try:
        opener = gzip.open if path.lower().endswith(".gz") else open
        with opener(path, "rb") as fh:
            head = fh.read(5)
    except (OSError, EOFError):
        return None
    if head == b"LOCUS":
        return "genbank"
    if head[:4] == b"ABIF":
        return "abif"
    if head[:1] == b"@":
        return "fastq"
    if head[:1] == b">":
        return "fasta"
    return None


def _read_fastq_text(path: str) -> Iterator[tuple[str, str, str]]:
    with _open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            if not h.strip():
                continue
            seq = fh.readline().strip()
            fh.readline()                       # '+'
            qual = fh.readline().strip()
            yield (h[1:].split() or ["unnamed"])[0], seq.upper(), qual


def _read_fasta_reads(path: str, default_q: int) -> Iterator[tuple[str, str, str]]:
    qc = chr(default_q + 33)
    for name, seq in read_fasta(path):
        if seq:
            yield name, seq, qc * len(seq)


_SEQ_CHARS = re.compile(r"[^ACGTUNRYKMSWBDHV]")


def _read_seq_text(path: str, default_q: int) -> Iterator[tuple[str, str, str]]:
    """Plain-text Sanger sequence: one read, named after the file."""
    with _open(path) as fh:
        text = fh.read()
    if text.lstrip().startswith(">"):
        yield from _read_fasta_reads(path, default_q)
        return
    if text.lstrip().startswith("@"):
        yield from _read_fastq_text(path)
        return
    letters = re.sub(r"[^A-Z]", "", text.upper())
    seq = _SEQ_CHARS.sub("", letters).replace("U", "T")
    # a text file that is not mostly nucleotides (notes, a report) is not a read
    if seq and len(seq) >= 0.9 * len(letters):
        yield _stem(path), seq, chr(default_q + 33) * len(seq)


def _read_genbank(path: str, default_q: int) -> Iterator[tuple[str, str, str]]:
    """Each GenBank record as a read: named by its LOCUS, sequence from ORIGIN."""
    name, seq, in_seq = None, [], False
    with _open(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                parts = line.split()
                name, seq, in_seq = (parts[1] if len(parts) > 1 else _stem(path)), [], False
            elif line.startswith("ORIGIN"):
                in_seq = True
            elif line.startswith("//"):
                s = _SEQ_CHARS.sub("", "".join(seq).upper()).replace("U", "T")
                if s:
                    yield name or _stem(path), s, chr(default_q + 33) * len(s)
                name, seq, in_seq = None, [], False
            elif in_seq:
                seq.append(re.sub(r"[^A-Za-z]", "", line))


# ------------------------------------------------------------------ ABIF --

_ABIF_HEADER = struct.Struct(">4sH4sihhiiii")      # magic, version, root dir entry
_ABIF_ENTRY = struct.Struct(">4sihhiiii")          # name, number, type, size, n, bytes, off, h


def parse_abif(data: bytes) -> dict[tuple[str, int], bytes]:
    """Raw ABIF directory: {(tag name, tag number): data bytes}.

    The format (Applied Biosystems, 2009 spec) is a 128-byte header whose
    root entry points at a directory of 28-byte entries; an entry's data sits
    at its offset, or inside the offset field itself when 4 bytes or fewer.
    """
    if len(data) < _ABIF_HEADER.size or data[:4] != b"ABIF":
        raise ValueError("not an ABIF file (missing 'ABIF' magic)")
    (_m, _ver, _name, _num, _etype, _esize, n_entries, _dsize,
     dir_off, _h) = _ABIF_HEADER.unpack_from(data, 0)
    tags: dict[tuple[str, int], bytes] = {}
    for i in range(n_entries):
        off = dir_off + i * _ABIF_ENTRY.size
        if off + _ABIF_ENTRY.size > len(data):
            break
        name, number, _t, _s, _n, dsize, doff, _hh = _ABIF_ENTRY.unpack_from(data, off)
        if dsize <= 4:
            raw = data[off + 20: off + 20 + dsize]
        else:
            raw = data[doff: doff + dsize]
        tags[(name.decode("latin-1"), number)] = raw
    return tags


def read_abif(path: str) -> tuple[str, str, str]:
    """(name, seq, qual) from an .ab1/.abi Sanger trace.

    Uses the edited basecalls (PBAS 2) and their qualities (PCON 2), falling
    back to the original calls (PBAS 1 / PCON 1). Without quality values the
    bases get FASTA_DEFAULT_Q. The read is named after the file, so a folder
    of clones gives one read per clone with a recognisable name.
    """
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rb") as fh:
        data = fh.read()
    try:
        tags = parse_abif(data)
    except (ValueError, struct.error):
        tags = {}
    seq_b = tags.get(("PBAS", 2)) or tags.get(("PBAS", 1))
    if seq_b is None:
        try:                                       # pragma: no cover - optional path
            from Bio import SeqIO
            rec = SeqIO.read(path, "abi")
            q = rec.letter_annotations.get("phred_quality")
            seq = str(rec.seq).upper()
            qual = "".join(chr(min(x, 93) + 33) for x in q) if q else \
                chr(FASTA_DEFAULT_Q + 33) * len(seq)
            return _stem(path), seq, qual
        except ImportError:
            pass
        raise ValueError(f"{path}: no basecalls (PBAS) in ABIF file")
    seq = seq_b.decode("latin-1").upper().replace("\x00", "")
    q_b = tags.get(("PCON", 2)) if ("PBAS", 2) in tags else tags.get(("PCON", 1))
    if q_b is None:
        q_b = tags.get(("PCON", 1)) or tags.get(("PCON", 2))
    if q_b is not None and len(q_b) >= len(seq):
        qual = "".join(chr(min(b, 93) + 33) for b in q_b[:len(seq)])
    else:
        qual = chr(FASTA_DEFAULT_Q + 33) * len(seq)
    return _stem(path), seq, qual


# ---------------------------------------------------------- paired reads --

def _phred_char(q: int) -> str:
    return chr(max(0, min(q, 93)) + 33)


def merge_pair(seq1: str, qual1: str, seq2: str, qual2: str, min_overlap: int = 20,
               max_mismatch_rate: float = 0.1, k: int = 12) -> tuple[str, str] | None:
    """Merge one overlapping read pair; None if they do not overlap.

    R2 is reverse-complemented and placed against R1 at the offset with the
    best overlap (at least `min_overlap` bases, at most `max_mismatch_rate`
    mismatches, scored matches - 4 x mismatches). Candidate offsets come from
    shared k-mers, so the cost is linear in read length. In the overlap, agreeing
    bases get the summed quality (capped at 60); disagreeing bases take the
    higher-quality call with quality |q1 - q2| (an 'N' at Q2 on a tie). When
    the fragment is shorter than the reads, the read-through beyond the other
    read's start (adapter) is dropped.
    """
    s2 = _revcomp(seq2.upper())
    q2 = [ord(c) - 33 for c in qual2[::-1]]
    s1 = seq1.upper()
    q1 = [ord(c) - 33 for c in qual1]
    n1, n2 = len(s1), len(s2)
    if n1 < min_overlap or n2 < min_overlap:
        return None
    k = min(k, min_overlap)
    index: dict[str, list[int]] = {}
    for i in range(n1 - k + 1):
        index.setdefault(s1[i:i + k], []).append(i)
    offsets = set()
    for j in range(0, n2 - k + 1, max(1, k // 2)):
        for i in index.get(s2[j:j + k], ()):
            offsets.add(i - j)
    best = None
    for o in offsets:
        a, b = max(0, o), min(n1, o + n2)
        L = b - a
        if L < min_overlap:
            continue
        mm = sum(1 for x in range(a, b) if s1[x] != s2[x - o])
        if mm > max_mismatch_rate * L:
            continue
        score = (L - mm) - 4 * mm
        if best is None or score > best[0]:
            best = (score, o)
    if best is None:
        return None
    o = best[1]
    a, b = max(0, o), min(n1, o + n2)
    seq, qual = [], []
    if o > 0:                                   # R1 sequence before R2 starts
        seq.append(s1[:o])
        qual.append(qual1[:o])
    for x in range(a, b):
        b1, b2, x1, x2 = s1[x], s2[x - o], q1[x], q2[x - o]
        if b1 == b2:
            seq.append(b1)
            qual.append(_phred_char(min(x1 + x2, 60)))
        elif x1 == x2:
            seq.append("N")
            qual.append(_phred_char(2))
        else:
            seq.append(b1 if x1 > x2 else b2)
            qual.append(_phred_char(max(2, abs(x1 - x2))))
    if o + n2 > n1:                             # R2 sequence after R1 ends
        seq.append(s2[b - o:])
        qual.append("".join(_phred_char(q) for q in q2[b - o:]))
    return "".join(seq), "".join(qual)


def merge_pairs(r1_path: str, r2_path: str, min_overlap: int = 20,
                max_mismatch_rate: float = 0.1,
                stats: dict | None = None) -> Iterator[tuple[str, str, str]]:
    """Merge overlapping Illumina paired-end reads into single reads.

    R1 and R2 are read in step (the usual bcl2fastq/BCL Convert output order).
    Pairs that do not overlap by `min_overlap` bases within the mismatch
    tolerance are dropped: an SHM amplicon is designed to be spanned, and an
    unmerged pair would count the overlap twice. Pass a dict as `stats` to
    receive pairs_in / pairs_merged / pairs_unmerged.
    """
    n_in = n_ok = 0
    r2 = read_fastq(r2_path)
    for name, s1, q1 in read_fastq(r1_path):
        try:
            _n2, s2, q2 = next(r2)
        except StopIteration:
            break
        n_in += 1
        m = merge_pair(s1, q1, s2, q2, min_overlap, max_mismatch_rate)
        if m is None:
            continue
        n_ok += 1
        yield re.sub(r"/[12]$", "", name), m[0], m[1]
    if stats is not None:
        stats.update(pairs_in=n_in, pairs_merged=n_ok, pairs_unmerged=n_in - n_ok)


_PAIR_RE = re.compile(r"^(.*?)([._-])(R?)([12])((?:[._-].*)?)$")


def find_pairs(paths: list[str]) -> tuple[list[tuple[str, str]], list[str]]:
    """Split file paths into ([(R1, R2), ...], singles).

    Two files pair when their names (minus extension) differ only in the read
    number of an _R1/_R2 (or _1/_2, also with . or -) token, e.g.
    sample_S1_L001_R1_001.fastq.gz and sample_S1_L001_R2_001.fastq.gz.
    """
    keyed: dict[tuple, dict[str, str]] = {}
    for p in paths:
        base = os.path.basename(p)
        stem, ext = base, ""
        low = base.lower()
        for e in sorted(READ_EXTENSIONS, key=len, reverse=True):
            for full in (e + ".gz", e):
                if low.endswith(full):
                    stem, ext = base[:-len(full)], full
                    break
            if ext:
                break
        m = _PAIR_RE.match(stem)
        if not m:
            continue
        key = (os.path.dirname(p), m.group(1), m.group(2), m.group(3).upper(),
               m.group(5), ext)
        keyed.setdefault(key, {})[m.group(4)] = p
    pairs, paired = [], set()
    for d in keyed.values():
        if "1" in d and "2" in d:
            pairs.append((d["1"], d["2"]))
            paired.update(d.values())
    singles = [p for p in paths if p not in paired]
    return sorted(pairs), singles


def split_pair_spec(path: str) -> tuple[str, str] | None:
    """'R1.fastq.gz+R2.fastq.gz' -> (R1, R2); None if `path` is not such a spec."""
    if PAIR_SEPARATOR not in path or os.path.exists(path):
        return None
    r1, _, r2 = path.partition(PAIR_SEPARATOR)
    r1, r2 = r1.strip(), r2.strip()
    if r1 and r2 and os.path.exists(r1) and os.path.exists(r2):
        return r1, r2
    return None


# ---------------------------------------------------------------- dispatch --

def _read_one(path: str, default_q: int, fmt: str | None = None):
    fmt = fmt or sniff_format(path)
    if fmt == "bam":
        yield from read_bam(path)
    elif fmt == "abif":
        yield read_abif(path)
    elif fmt == "fasta":
        yield from _read_fasta_reads(path, default_q)
    elif fmt == "seq":
        yield from _read_seq_text(path, default_q)
    elif fmt == "genbank":
        yield from _read_genbank(path, default_q)
    elif fmt == "fastq":
        yield from _read_fastq_text(path)
    else:
        raise ValueError(f"{path}: unrecognised read file (expected one of "
                         f"{', '.join(READ_EXTENSIONS)}, optionally .gz)")


def _read_dir(path: str, default_q: int):
    files = []
    for name in sorted(os.listdir(path)):
        full = os.path.join(path, name)
        if name.startswith(".") or not os.path.isfile(full):
            continue
        files.append(full)
    pairs, singles = find_pairs([f for f in files if sniff_format(f) == "fastq"])
    paired = {p for pr in pairs for p in pr}
    for r1, r2 in pairs:
        yield from merge_pairs(r1, r2)
    for f in files:
        if f in paired:
            continue
        fmt = sniff_format(f)
        if fmt is None:
            continue                            # notes, PDFs, ... in a clone folder
        yield from _read_one(f, default_q, fmt)


def read_fastq(path: str, default_q: int = FASTA_DEFAULT_Q) -> Iterator[tuple[str, str, str]]:
    """Every read of one sample as (name, seq, qual), whatever the input format.

    See the module docstring for what `path` may be. `default_q` is the Phred
    given to bases of inputs without qualities (FASTA, .seq, GenBank).
    """
    path = os.fspath(path)
    if os.path.isdir(path):
        yield from _read_dir(path, default_q)
        return
    pair = split_pair_spec(path)
    if pair:
        yield from merge_pairs(*pair)
        return
    yield from _read_one(path, default_q)


def phred(qual: str) -> list[int]:
    return [ord(c) - 33 for c in qual]


def expected_error_rate(qual: str) -> float:
    """Mean per-base error probability = sum(10^(-Q/10)) / length (Q20-eq at 0.01)."""
    if not qual:
        return 1.0
    return sum(10 ** (-q / 10.0) for q in phred(qual)) / len(qual)


def mott_trim(qual: str, cutoff: float = 0.01) -> tuple[int, int]:
    """Modified-Mott quality trim: [start, end) of the best-scoring segment.

    Each base scores cutoff - 10^(-Q/10) (positive above the quality the cutoff
    corresponds to: 0.01 = Q20, 0.05 = Q13); the kept segment is the contiguous
    run with the highest total. This is phred/Biopython's trace trimming and
    removes the noisy start and decaying end of a Sanger read in one pass.
    Returns (0, 0) when no base clears the cutoff.
    """
    best_sum = run = 0.0
    best = (0, 0)
    start = 0
    for i, c in enumerate(qual):
        s = cutoff - 10 ** (-(ord(c) - 33) / 10.0)
        if run <= 0:
            run, start = s, i
        else:
            run += s
        if run > best_sum:
            best_sum, best = run, (start, i + 1)
    return best

# Paths for regenerating the manuscript figures with the corrected Sanger data.
# Nanopore inputs are read, unchanged, from the original analysis folder; the
# corrected Sanger calls come from ../sanger_calls (built from the organised .ab1
# traces); figures are written to ../figures. Nothing in the original folders is
# modified.
import os
BASE        = ("/Users/abdulmateenaderinto/Documents/Projects/2026/PhD - Fundamental Immunology/"
               "Martin Lab/Paper #1 - Nanopore sequencing/June 2nd analyses")
DATA_JH4    = os.path.join(BASE, "July 13th JH4")
DATA_SANGER = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sanger_calls")
OUT_DIR     = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
REF_FASTA   = os.path.join(DATA_JH4, "DHJGKB_2_Wildtype_2_C57BL6_sub1000 assembled to Jh4.fasta")
WINDOW      = (101, 629)
FULL_REGION = (90, 1230)
DEPTH       = "sub500"
os.makedirs(OUT_DIR, exist_ok=True)

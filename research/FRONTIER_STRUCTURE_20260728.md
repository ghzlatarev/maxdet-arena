# Order-23 frontier structure and local H-equivalence audit

**Audit date:** 2026-07-28

This note separates three claims:

1. exact arithmetic implied by the arena contract and published bounds;
2. structure read from the published order-23 Gram matrix;
3. a reproducible classification of the **local** matrix corpus.

The six-class result below is not a classification of all order-23 record
matrices. It says only what is represented by the listed local artifacts.
A continued-search addendum records a seventh local class found after that
input snapshot was frozen.

## Exact score lattice

The published comparison matrix verifies at

```text
D0 = 2,779,447,296,000,000
   = 2^22 * 662,671,875
   = 2^22 * 3 * 5^6 * 67 * 211.
```

Every order-23 sign determinant is divisible by `2^22 = 4,194,304`. Write
`|det A| = 2^22 d`. The first score lattice point strictly above the published
record is therefore

```text
d = 662,671,876
|det A| = 2,779,447,300,194,304.
```

The order-specific Ehlich upper bound used by the verifier is

```text
U = 2^22 * 3 * 5^6 * 675 * sqrt(505)
U^2 = 8,894,085,385,420,800,000,000,000,000,000.
```

Set `C = 3 * 5^6 * 675 = 31,640,625`. Combining this bound with determinant
divisibility gives the exact integer-lattice ceiling

```text
d <= floor(C * sqrt(505)) = 711,034,613
|det A| <= 2,982,295,321,444,352.
```

Thus any strict improvement must satisfy

```text
662,671,876 <= d <= 711,034,613.
```

These are necessary score constraints, not a list of attainable determinants.
A candidate Gram matrix must additionally be positive definite, have diagonal
23, obey the sign-matrix congruences, have square determinant
`(2^22 d)^2`, and admit an exact factor in `{-1,+1}^{23x23}`.

The original 2003 construction, Orrick's 2006 enumeration paper, the 2021
survey, and the archived order-23 table all retain `D0` as the published
comparison value. No stronger published construction was found in the sources
surveyed through this audit date. That literature result is not a claim about
unpublished computations.

## The useful Gram shell

After row switching, the published row Gram can be written

```text
G = 24 I - J + 4 B,
```

where `B` is the adjacency matrix of the pairs with inner product `+3`;
nonedges have inner product `-1`.

For `n = 23`, the ideal Ehlich block partition is

```text
3 + 4 + 4 + 4 + 4 + 4,
```

so its defect graph is `K3 disjoint-union 5 K4`, with 33 edges. Its determinant
is exactly the squared Ehlich expression,

```text
det(G_ideal)
  = 505 * (2^22 * 3 * 5^6 * 675)^2
  = 2^44 * 3^8 * 5^17 * 101.
```

This is not a square, so `G_ideal` cannot equal `A A^T` for a nonsingular
integer matrix `A`.

The published record Gram is an especially structured 12-edge augmentation of
that ideal graph. In one 1-based relabeling, take the base cliques

```text
{1,2,3}
{4,5,6,7} {8,9,10,11} {12,13,14,15}
{16,17,18,19} {20,21,22,23}.
```

Add the three complete bipartite connectors

```text
{1,3} x {6,7}
{2,3} x {10,11}
{1,2} x {14,15}.
```

Each connector contributes four edges. Two `K4` blocks remain untouched. The
result has 45 edges and degree multiset

```text
3^14, 5^6, 6^3.
```

Its graph automorphism group has order `442,368 = 2^14 * 3^3`. This identifies
a high-value exact search target: enumerate 12-edge augmentations and adjacent
shells of `K3 disjoint-union 5 K4` modulo automorphisms. The earlier local
`ideal + at most 4 edges` screen did not reach the shell containing the
published record.

## Local H-equivalence audit

### Method

For every pivot `(r,c)`, dephase a sign matrix `A` by

```text
N[i,j] = A[i,j] A[r,j] A[i,c] A[r,c].
```

Pivot row `r` and pivot column `c` are then positive. Delete them and encode
the `-1` entries of the remaining `22 x 22` matrix as a bipartite graph with
separate row and column colors. A color-preserving nauty certificate removes
the remaining row and column permutations. Taking the lexicographically
minimum certificate over all `23^2 = 529` pivots gives a complete certificate
for signed row/column permutation (H-) equivalence.

The audit also repeats the calculation after transposition and groups by
`min(cert(A), cert(A^T))`. Certificate hashes below use
`pynauty==2.8.8.1`.

### Inputs

The frozen local-corpus input set in `h_equivalence_audit.py` contains these
named artifacts plus every
`runs/direct-search/gram-factors/*/*.matrix.txt`:

- the published reference and the recovered pre-April-2003 factor;
- H24-derived class-14, class-9, and class-51 frontier factors;
- three depth-12 beam ties;
- four shell-MILP factors (seeds 29761, 29762, 29768, and 29773);
- sphere ties `6e5f51f4...`, `eb29a461...`, and `85a63a49...`;
- 40 Gram-factor files, including 14 byte-for-byte duplicates.

Every unique input is independently checked to have exact absolute determinant
`D0` before classification.

### Result

The 55 paths contain 41 byte-distinct frontier matrices and 14 duplicate
files. They form exactly six H-classes:

| H-certificate SHA-256 | Count | Local members |
|---|---:|---|
| `133097154bd369ae928c5c712f50745afd75f1d4265e2137dc57475cdb0c6d99` | 4 | H24 class 9; third depth-12 tie; MILP 29761; sphere tie `eb29a461...` |
| `5f3d7a03b4434c9e4432b084175480465bcbd410d8fba2fc5d360dd94667a63f` | 2 | MILP 29768 and 29773 |
| `6114b83b5af6eadbabe64c941e9c6706293aec09490c33d6407e78559d87e482` | 5 | H24 class 14; first and second depth-12 ties; sphere ties `6e5f51f4...` and `85a63a49...` |
| `9035bdf2a85b8a2a600a76c6d55af36f627327694e904e7b483e88993716c91b` | 1 | H24 class 51 |
| `b584c923ea12af6634bb5681f1eb903196e5a7cedadf5d1ea321bc63613c986e` | 1 | MILP 29762 |
| `db2cddf4b8f12da99a32b1563689ad693b89d0a9ae343a4aedde56361fdc5a81` | 28 | reference; pre-April factor; 26 byte-distinct randomized Gram factors |

Allowing transposition still gives six classes; no two of the six local
H-classes merge. All 41 matrices have the same normalized row-Gram defect
graph:

```text
Gram certificate SHA-256:
e641e83822fae412ec3e1e4d8aa355f015d204f10fa98adde97c5da508fbf554

edges: 45
automorphism-group order: 442,368
```

The factor-generation campaign has therefore found several genuinely distinct
factor classes, but it has not diversified beyond one Gram equivalence class.

Orrick's 2006 Table 1 reports at least 14 inequivalent order-23 matrices at
this determinant. Consequently, this six-class local corpus omits at least
eight reported H-classes. Recovering those factors, or exhaustively enumerating
factors of the known Gram modulo its large automorphism group, should precede
more search time from raw-hash-diverse copies of already represented classes.

### Continued-search addendum: a seventh local H-class

The exact-pair QUBO trust-region pilot retained a Hamming-12 frontier tie:

```text
path:
runs/qubo-trust-pilot-20260728-seed31003/best-proposal.matrix.txt

absolute determinant:
2,779,447,296,000,000

raw SHA-256:
b386a8714a7bba8e2e88796f136e607ae56402f45695b31bc52f87b57ef0ca52

H and HT certificate SHA-256:
b64c33090c9aac4dc7386213917cec66ac28d2da9685319b585067e95d0d63f6
```

Adding this matrix to the frozen audit changes the local summary from 41 to
42 byte-distinct matrices and from six to seven H-classes. Allowing
transposition also gives seven classes. The normalized row-Gram certificate
does not change, so all 42 factors remain in the same Gram class.

This is a new class relative to the frozen local files, not a literature
novelty claim and not a score improvement. Orrick's reported set is still
larger than the seven classes represented locally. The concise audit artifact
is
`runs/qubo-trust-pilot-20260728-seed31003/h-equivalence-audit.json`.

## Reproduce

The classifier is outside the verifier and does not modify matrices:

```sh
python3 -m venv /tmp/maxdet-h-audit
/tmp/maxdet-h-audit/bin/pip install pynauty==2.8.8.1
/tmp/maxdet-h-audit/bin/python research/h_equivalence_audit.py \
  --local-corpus > /tmp/maxdet-h-audit.json
```

Expected summary:

```text
input_file_count      55
duplicate_file_count  14
unique_matrix_count   41
h_class_count          6
ht_class_count         6
gram_class_count       1
```

To reproduce the continued-search addendum, append the candidate explicitly:

```sh
/tmp/maxdet-h-audit/bin/python research/h_equivalence_audit.py \
  --local-corpus \
  runs/qubo-trust-pilot-20260728-seed31003/best-proposal.matrix.txt
```

The resulting counts are 56 input paths, 42 byte-distinct matrices, seven
H-classes, seven classes allowing transposition, and one normalized row-Gram
class.

The `runs/` artifacts are local research evidence and are not guaranteed to
ship in a clean public clone. The script's explicit path manifest and JSON
output make the classification reproducible when that corpus is available.

## Sources

- W. P. Orrick, B. Solomon, D. Dowdeswell, and W. Smith,
  [New lower bounds for the maximal determinant problem](https://arxiv.org/abs/math/0304410)
  (2003): record matrix and Gram matrix.
- W. P. Orrick,
  [On the enumeration of some D-optimal designs](https://arxiv.org/abs/math/0511141)
  (2006): Table 1 reports at least 14 inequivalent order-23 record matrices.
- P. Browne, R. Egan, F. Hegarty, and P. Ó Catháin,
  [A survey of the Hadamard maximal determinant problem](https://arxiv.org/abs/2104.06756)
  (2021): Ehlich block bound and historical determinant table.
- Orrick's archived
  [order-23 determinant page](https://web.archive.org/web/20200219170713id_/http://www.indiana.edu/~maxdet/d23.html):
  normalized score, printed Gram, factor, and comparison ratio.
- R. P. Brent,
  [Finding D-optimal designs by randomised decomposition and switching](https://arxiv.org/abs/1112.4671):
  exact Gram-factor search context.

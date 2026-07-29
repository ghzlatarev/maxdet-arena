# Fast exact principal-minor cube search, 2026-07-28

## Result

The new evaluator exhausts an arbitrary cube of up to 27 distinct entry
flips.  A 200-cube campaign checked exactly:

```text
200 × 2^27 = 26,843,545,600 matrices
best |det| = 2,779,447,296,000,000
strict gains = 0
```

The verified frontier is unchanged.  The batch recovered 28 nonzero frontier
tie assignments.  A bounded rerun captured all 28, deduplicated them to 10 raw
matrices, and independently passed Bareiss and `./arena verify` on every one.
All 10 were already present in the known 24-state neutral union.

## Exact reduction

For an invertible base matrix \(A\) and distinct flip coordinates
\((r_i,c_i)\), set \(\delta_i=-2A_{r_i,c_i}\).  The determinant lemma gives

\[
\frac{\det(A_S)}{\det(A)}
  = \det(M[S,S]),\qquad
M_{ij}=[i=j]+\delta_i(A^{-1})_{c_i,r_j}.
\]

All \(2^m\) candidates are therefore the principal minors of one \(m\times m\)
matrix.  `fast_principal_cube.cpp` computes them by a breadth-first Schur
recurrence over the field \(p=2^{32}-5\), using one batch inversion per level.
An exact zero pivot is temporarily shifted by one and corrected in reverse
dependency order.  The final scalar level needs no inversion.

For every order-23 sign matrix, the determinant is divisible by \(2^{22}\),
and

```text
floor(sqrt(23^23) / 2^22) = 1,089,457,290 < p/2.
```

The centered residue of `det(A)/2^22` consequently recovers the unique signed
integer determinant inside the Hadamard bound.  Every strict promotion, final
best, and emitted tie artifact crosses an independent integer-Bareiss check.
Matrix, log, and report files use fsync plus atomic rename.

The recurrence is based on Griffin and Tsatsomeros, “Principal minors, Part I:
A method for computing all the principal minors of a matrix,” *Linear Algebra
and its Applications* 419 (2006),
<https://doi.org/10.1016/j.laa.2006.04.008>.  This implementation was written
for the arena; no archived third-party source was copied.

## Validation

The strict build gate is:

```sh
c++ -std=c++20 -O3 -march=native \
  -Wall -Wextra -Wshadow -Werror -pedantic \
  research/fast_principal_cube.cpp \
  -o build/research/fast_principal_cube
```

`--self-test` covers:

- 256 random and deliberately singular principal-minor systems through
  dimension eight, compared with direct modular determinants;
- 72 random end-to-end sign-matrix entry-flip cubes through dimension eight,
  with every assignment compared against exact Bareiss;
- explicit singular candidate and zero-pivot cases;
- duplicate-coordinate rejection;
- centered quotient recovery at both Hadamard-bound endpoints;
- a SHA-256 known-answer test.

The same suite passes ASan and UBSan.  The repository's 39 trusted tests pass,
and the retained 27-bit result passes `./arena verify`.

An independent scratch audit additionally covered all 19,683 three-by-three
matrices over `{0,1,2}`, 17,600 random systems through order eight, and 2,500
random end-to-end flip cubes.  It confirmed the BFS mask ordering, the
determinant-lemma orientation, and the necessity of reverse zero-pivot
correction.

## Controlled benchmark

The legacy exact run between two 24-entry-distant frontier ties used one
Bareiss determinant per state:

```text
assignments: 16,777,216
elapsed:     356.410855 s
ties:        3
```

The new evaluator reproduced the same three ties:

```text
assignments: 16,777,216
elapsed:     0.241662 s
throughput:  69,424,303 candidates/s
max RSS:     378,175,488 bytes
speedup:     1,474.8×
```

Evidence:

```text
runs/direct-search/fast-principal-cube/known24-20260728/report.json
```

## Curated 27-bit cube

The first new support combined the verified class-9 `A0` neutral switch with a
15-entry row/column-disjoint transverse matching.  The transverse coordinates
were ranked by their best exact one-flip score at either neutral endpoint and
excluded the known 72-coordinate six-generator support.

```text
assignments: 134,217,728
elapsed:     2.423214 s
max RSS:     1,786,200,064 bytes
frontier ties: 2
strict gains: 0
```

Evidence:

```text
research/class9_a0_transverse27.coords.txt
runs/direct-search/fast-principal-cube/class9-a0-transverse27-20260728/report.json
```

## Mixed 200-cube campaign

`fast_cube_batch.py` generated exactly 40 affine-distinct supports from each
of five strategies:

- top exact one-flip scores;
- strongest exact pair score and pair-rescue linkage;
- recovered neutral switches plus transverse coordinates;
- random bipartite supports with all 23 rows and columns represented and line
  degree at most two;
- overlapping destroy/repair supports built from linked alternatives.

The four neutral-network starts cover both established local H-classes.  The
fifth start is the independently verified QUBO tie with raw SHA-256
`b386a871...`, classified separately as the newly observed seventh local
H/HT class.  Affine-cube fingerprints encode fixed signs and free coordinates,
so choosing a different endpoint inside the same cube cannot evade
deduplication.

```text
supports:                    200
unique affine fingerprints: 200
H0 / H1 / H2 cubes:          80 / 80 / 40
assignments:                 26,843,545,600
engine seconds:              383.359945
wall seconds:                390.922356
wall throughput:             68,667,205 candidates/s
maximum RSS:                 2,037,907,456 bytes
zero-pivot corrections:      0
promotions:                  0
nonzero frontier ties:       28
```

The real batch did not need a zero-pivot correction; that path is exercised by
the differential tests.  The 28 tie masks occurred in 26 neutral/transverse or
destroy/repair cubes, all on H0/H1.  No H2 cube tied anywhere except its start.

Evidence:

```text
runs/direct-search/fast-principal-cube/batch200-mixed-h012-20260728/manifest.json
runs/direct-search/fast-principal-cube/batch200-mixed-h012-20260728/aggregate-report.json
runs/direct-search/fast-principal-cube/batch200-mixed-h012-20260728/tie-audit/report.json
```

## Artifact hashes

```text
8280a872d3197caa08c2d431b5cdc04f95514acc5af787755ee2d2efb329563c  research/fast_principal_cube.cpp
118efc94299d9220d0cbabf9e5181c93885cc474ea932c6cc95c2020217af975  research/fast_cube_batch.py
a0904e62d902dbe315fb43d207e1e86dc1aff8454f85ed0c9ff0940894a7111f  research/fast_cube_tie_audit.py
a31b7c87e8a234fe239850a1c2c69c2687af8550c1eafcbac7f581471386cadf  research/class9_a0_transverse27.coords.txt
c943789a1fcac5df585c2a912ed7cf027a2ed41cfe6ee21eeb79629475037411  known24-20260728/report.json
ab310823916b8b4560922b4184166d3e115a4cc6ef5e65265b3bd1d569c1e67c  class9-a0-transverse27-20260728/report.json
176e37e3575f1c64e8cc6fcad4344689edfb261411d30f40753b332642b755b1  batch200-mixed-h012-20260728/manifest.json
2ad60931d99d5bd5d5874b7969dcaa46ee05f87e088e95ad1bea6256a52f75e3  batch200-mixed-h012-20260728/aggregate-report.json
ff86b13e910ed1f3765ab47a12473eeaa4b71ad9a209476d1bff7c1f177c1946  batch200-mixed-h012-20260728/tie-audit/report.json
```

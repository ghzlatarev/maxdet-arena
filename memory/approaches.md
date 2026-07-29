# Shared research memory

Keep this short and evidence-based. Submission-specific details belong in each
submission’s `notes.md`.

## Published comparison point

The order-23 matrix printed by Orrick, Solomon, Dowdeswell, and Smith (2003)
reproduces exactly:

```text
|det| = 2,779,447,296,000,000
sign-normalized sha256 = c0ea58d361945b20dad78bddb3fd93c0810b762e527e8781e87d6db5e86d993a
```

The repository matrix is a byte-for-byte match after mechanically translating
the `+`/`-` rows in the authors’ archived `matData.tex`. That source reports
the tighter order-23 Ehlich upper bound
`2^22 * 3 * 5^6 * 675 * sqrt(505)` and a reference ratio of `0.931983`.
The 2021 survey’s `0.7091` uses the generic Barba bound and explicitly notes
that its table does not use the best known bound for orders at most 59.

## Prelaunch dogfooding

- The deterministic `J - 2I` genesis matrix verifies at `88,080,384`.
- A dependency-free random first-improvement search reached
  `2,088,024,410,161,152` in ten seconds on the initial development machine.
  This is useful as a smoke test but remains below the published comparison
  point.
- Direct execution initially failed to find the local `maxdet` package. The
  starter now pins its repository import root and is tested from a clean clone.
- The first strict parser rejected a harmless extra trailing newline. The
  contract now accepts surrounding ASCII whitespace while rejecting Unicode
  whitespace and every non-domain character.
- Raw hashes differ when semantically identical whitespace changes. Receipts
  therefore expose both a raw file hash and a deterministic sign-normalized
  matrix hash.
- Long native searches initially emitted nothing between start and finish.
  Search logs now include periodic heartbeats, elapsed time, accepted moves,
  seed, mode, and exact incumbent score.
- Submission preparation now validates score and metadata before publishing an
  atomic bundle. A tie, regression, bad slug, or blank method leaves no partial
  directory behind.
- CI derives the effective frontier from the published floor, exactly
  reproduced arena checkpoint, and every verified trusted submission rather
  than relying only on a manually updated display number.
- Exact-score cross-checks now cover both the matrix and Gram matrix modulo
  three primes; their combined modulus uniquely identifies a determinant
  inside the Hadamard bound.
- Whole-row/column coordinate moves use inverse cofactors to propose the
  determinant-maximizing sign pattern for one line. Early inverse-only gain
  tests entered false-progress loops near numerical tolerance, so every block
  move now proves a strict exact determinant increase before acceptance.
- An exhaustive exact audit of all 24,673,089 perturbations at Hamming radius
  one, two, or three from the published matrix found no strict improvement
  (555.77 seconds in the prelaunch run). This establishes only radius-three
  local optimality for that representation, not global optimality.
- A completion-counted exact two-line audit jointly optimized all 253 row
  pairs and all 253 column pairs. It evaluated 2,122,317,824 sign assignments
  in 61.735 seconds and found no strict improvement. This establishes only
  local optimality under replacement of two complete parallel lines.
- Time-bounded finite searches must not be called exhaustive without a
  completion count. An early pair audit reached only 471 of 506 pairs, so the
  research tool gained a pass-counted mode before the retained audit.

## Search directions

- Rank-one determinant updates make single-entry local search cheap, but inverse
  drift must be periodically rebuilt and every candidate must be exactly
  reverified.
- Try multiple basins: random restarts, perturbations of the 2003 matrix,
  annealing, tabu search, row-pattern replacements, Gram-matrix methods, and
  SAT/constraint formulations.
- For a fixed set of 22 rows, the determinant is linear in the remaining row.
  The best sign row is the sign pattern of its cofactor vector; the same holds
  for columns. Alternating exact-accepted line replacements is a useful search
  primitive and a stronger local-optimality test than single-entry flips.
- For two fixed row positions, second cofactors reduce joint replacement to
  `max_x sum_j |sum_i x_i C_ij|`. Fixing one redundant sign leaves `2^22`
  assignments, which can be traversed cheaply in Gray-code order. Kicked
  exploration must keep its lower-scoring working state separate from the
  atomically checkpointed global incumbent.
- The published matrix has 208 row pairs with inner-product magnitude `1` and
  45 with magnitude `3`; no larger off-diagonal Gram magnitude occurs. The
  graph on the 45 magnitude-3 pairs has degree distribution
  `3^14, 5^6, 6^3`. Structured moves that preserve or deliberately rewire this
  sparse Gram-defect graph may be more promising than blind entry flips.
- A retained 12-entry-kick/two-line basin reached exact determinant
  `2,709,848,064,000,000 = 2^22 * 3 * 5^6 * 7 * 11 * 179`, or `90.8645%` of
  the order-23 Ehlich bound. Its row Gram matrix has 206 off-diagonal pairs of
  magnitude `1` and 47 of magnitude `3`. It is below the published comparison
  point but is a useful distinct, reproducible private checkpoint.
- Replaying the retained 24-entry basin improved the private checkpoint to
  `2,726,756,352,000,000 = 2^22 * 3^3 * 5^6 * 23 * 67`, or `91.4315%` of the
  Ehlich bound. Its row Gram matrix has 202 off-diagonal pairs of magnitude
  `1` and 51 of magnitude `3`; its exact receipt begins `0922f725`.
- Seed 2403 independently reached the same score with a different
  sign-normalized hash. Enumerating pivot normalizations and checking the
  resulting 22+22 vertex bipartite graphs found a row/column
  sign-and-permutation equivalence (pivot row 2, column 1 in one-based
  indexing). Do not count these as separate constructions.
- All three 7,200-second pair-kick seeds reached the same exact
  `2,726,756,352,000,000` score in a perturbed basin. A deterministic
  900-second seed-2402 replay retained it. Exact line ascent from that replay
  reached the published reference score in 0.151 seconds, so treat this state
  as a useful high waypoint in the reference basin rather than a separate
  line-local optimum.
- The completed private campaign used 32.42 single-core hours across
  first-improvement, hill, annealing, hybrid, coordinate, exact block,
  radius-three audit, and exact two-line search. None beat the published
  comparison point. The strongest purely random-coordinate run reached
  `2,391,226,070,335,488`, far below the kicked exact-pair basin.
- Sign flips and row/column permutations create large duplicate families. Do
  not count them as progress.

## Frontier-portal campaign (2026-07-29)

- Polar/Douglas--Rachford projection with the Ehlich target spectrum found two
  frontier H/HT classes absent from the frozen local corpus. The explicit
  local union grew from 8 to 10 H/HT classes (14 H classes) but remained one
  normalized row-Gram class. This is local-corpus novelty only.
- The substantive projection run exact-scored 224,740 distinct sign matrices,
  retained 155 search-discovered frontier ties (163 including inputs), and
  found no strict win. Flat-spectrum guidance found no new local class; the
  Ehlich geometry was the productive representation.
- Exact entry-flip radius at most three is closed around six previously
  uncovered portals and both new portals: 197,384,712 center/subset
  assignments total, with zero frontier ties and zero strict wins.
- Two exact integer functionals reject 14 of the 20 fixed-Gram small-shell
  triples. The six surviving orientations are the observed `(0,2,4)` orbit;
  this theorem does not classify all factorizations within that orbit.
- Exact pairwise alignment under `Aut(G)` gives Hamming distances
  `{32,46,47}` among the ten local H/HT portals. Four structurally different
  distance-32 affine connector cubes were exhaustively enumerated
  (`4 * 2^32 = 17,179,869,184` assignments); each contained only its two
  frontier endpoints, with no interior tie or strict win. Direct portal-to-
  portal interpolation is therefore a valley in these pinned geometries, not
  an observed above-frontier ridge.

# Experimental research tools

Nothing in this directory is part of the trusted submission verifier. Agents
may replace it freely; only `candidate/matrix.txt` crosses the boundary.

## Native search

Build the dependency-free C++ research binary:

```sh
mkdir -p build/research
c++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  research/fast_search.cpp -o build/research/fast_search
```

Example:

```sh
build/research/fast_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/block.matrix.txt \
  --log runs/block.jsonl \
  --mode block \
  --seed 909 \
  --seconds 3600 \
  --heartbeat-seconds 60

./arena verify runs/block.matrix.txt
```

Modes:

- `hill`: greedy one-entry ascent with random restarts;
- `anneal`: inverse-guided single-entry annealing;
- `hybrid`: greedy ascent plus small kicks;
- `coordinate`: random-restart whole-row/column coordinate ascent;
- `block`: kicks from the incumbent followed by exact-accepted line moves;
- `audit2` / `audit3`: exhaustive exact Hamming-neighborhood audits.

Floating-point inverses propose moves only. Candidate promotions and every
whole-line acceptance use integer Bareiss determinants; `./arena verify` remains
the authority and adds Gram, modular, bounds, divisibility, and hash checks.

Compiled binaries and `runs/` logs are intentionally ignored. Preserve durable
results as a small record or submission artifact, not as an executable.

## Valley-crossing searches

`reactive_tabu.cpp` ranks all 529 single-entry flips with inverse updates and
takes the best admissible move, including downhill moves. `symmetric_tabu.cpp`
instead preserves `A = A^T` and also proposes disjoint symmetric edge pairs.
The latter accepts `--restart-iterations N` (`0` disables restarts).

`triple_line_search.cpp` perturbs two rows or columns by one to four entries,
then replaces a third line with the determinant-maximizing cofactor-sign
pattern. It alternates orientations, accepts a controlled number of downhill
moves, and periodically kicks from its exact incumbent.

Build any of them with the same warning gate:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/triple_line_search.cpp -o build/research/triple_line_search

build/research/triple_line_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/triple.matrix.txt \
  --log runs/triple.jsonl \
  --seed 28901 --seconds 3600 --heartbeat-seconds 60
```

All floating-point work is proposal-only. Bareiss integers govern acceptance
and atomic best-checkpoint promotion.

## Fixed-radius and multiflip searches

`hamming_sphere_tabu.cpp` keeps exactly `r` entries flipped from a supplied
matrix and walks by one-out/one-in exchanges. Rank-two inverse updates score a
full or sampled exchange neighborhood; accepted states are checked by exact
Bareiss determinants. The tool supports reactive reverse-move tabu, controlled
downhill moves, deterministic restarts, multiple radii, an independent first
tie artifact, atomic checkpoints, and clean signal handling:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Wshadow -Werror -pedantic \
  research/hamming_sphere_tabu.cpp -o build/research/hamming_sphere_tabu

build/research/hamming_sphere_tabu \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/sphere-best.matrix.txt \
  --tie-output runs/sphere-tie.matrix.txt \
  --checkpoint runs/sphere-checkpoint.json \
  --log runs/sphere.jsonl \
  --score-floor 2779447296000000 \
  --radii 24,32,40 --seed 29901 --seconds 900 \
  --restart-iterations 1536 --swap-samples 0
```

`multiflip_beam.cpp` instead constructs bounded-cardinality entry sets one
flip at a time, retains a determinant-ranked beam, and applies randomized
fixed-cardinality refinements. `--max-per-row` and `--max-per-column` prevent
one line from consuming the move; depths through 24 are supported. Use
`--tie-output` to preserve the first non-baseline matrix at the score floor
without replacing the global checkpoint:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Wshadow -Werror -pedantic \
  research/multiflip_beam.cpp -o build/research/multiflip_beam

build/research/multiflip_beam \
  --start runs/frontier.matrix.txt \
  --output runs/multiflip-best.matrix.txt \
  --tie-output runs/multiflip-tie.matrix.txt \
  --log runs/multiflip.jsonl \
  --seed 29910 --runs 1000 --beam-width 12000 \
  --min-flips 12 --max-flips 12 \
  --max-per-row 5 --max-per-column 5 \
  --refine-states 192 --random-refinements 96 --swap-passes 3 \
  --seconds 900 --heartbeat-seconds 30
```

Both methods are stochastic searches. A completed wall-clock run is not an
exhaustive Hamming-ball statement, and every retained artifact still requires
`./arena verify`.

## Neutral-switch analysis

`neutral_cycle.py` infers six structured 12-entry generators from the
class-9 frontier ties and checks their full 64-state compound cube with exact
Bareiss determinants. `neutral_extension.py` adds an outside sphere-derived
pair and checks the resulting 256 selector states. `neutral_completion.py`
derives the complementary generators around that outside pair, tests every
orientation coupling, and emits the second exact neutral cycle.
`neutral_union.py` exhausts the full 13-generator union of both cycles and
their bridge. The analyzers write fresh reports and all frontier matrices
atomically:

```sh
python3 research/neutral_cycle.py \
  --base runs/class9.matrix.txt \
  --a0-endpoint runs/tie-a0.matrix.txt \
  --b0-endpoint runs/tie-b0.matrix.txt \
  --a1-endpoint runs/tie-a1.matrix.txt \
  --b1-endpoint runs/tie-b1.matrix.txt \
  --output-dir runs/neutral-cycle
```

`fixed_support_enum.cpp` then exhausts every labeled 12-edge mask with the
observed `(2,2,2,1,1,1,1,1,1)` row and column degrees on a fixed 9-by-9
support. Its determinant recovery is exact, using the universal `2^22`
divisibility, the Hadamard bound, and arithmetic modulo `2^32 - 5`; there is
no floating acceptance gate. The fixed-support result does not enumerate
other support choices or prove that the larger neutral network is complete.

## Exact recombinant cubes and crossover tabu

`crossover_tabu.cpp` restricts all changes to entries where two endpoint
matrices differ. `--mode exact` Gray-code enumerates every recombinant for a
small Hamming dimension; `--mode tabu` provides a persistent valley search
for larger endpoint differences:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Wshadow -Werror -pedantic \
  research/crossover_tabu.cpp -o build/research/crossover_tabu

build/research/crossover_tabu \
  --first runs/endpoint-a.matrix.txt \
  --second runs/endpoint-b.matrix.txt \
  --output runs/crossover-best.matrix.txt \
  --log runs/crossover.jsonl \
  --mode exact --seed 29916 --seconds 1200
```

Exact mode records complete assignment counts. Tabu mode uses inverse
proposals, exact aspirations and promotions, endpoint/random restarts, and
atomic checkpoints; it remains a stochastic negative search.

`fast_principal_cube.cpp` exhausts the same kind of entry-selector cube by
reducing every candidate determinant to one principal minor of a small
finite-field matrix.  For dimensions through 27 this replaces one Bareiss
elimination per assignment with a fast principal-minor recurrence.  Modular
quotient recovery is exact for order 23, zero pivots are corrected in reverse
dependency order, duplicate coordinates are rejected, and each promotion is
independently checked by integer Bareiss before an atomic checkpoint:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Wshadow -Werror -pedantic \
  research/fast_principal_cube.cpp \
  -o build/research/fast_principal_cube

build/research/fast_principal_cube --self-test

build/research/fast_principal_cube \
  --start runs/frontier.matrix.txt \
  --coordinates research/support.coords.txt \
  --output runs/cube/best.matrix.txt \
  --tie-output runs/cube/first-frontier-tie.matrix.txt \
  --log runs/cube/search.jsonl \
  --report runs/cube/report.json
```

Coordinate files contain distinct one-based `row column` pairs.  `--second`
may replace `--coordinates`; its differing entries become the cube support.
`--top-k N` retains the best nonzero strictly subfrontier masks in a
deterministic bounded heap; `--top-k-per-weight N` does the same separately
for every Hamming weight.  Optional `--top-k-output-dir` artifacts are
Bareiss-verified.  With all top-K options omitted, output and search behavior
remain unchanged.
The built-in test compares random and singular systems through dimension eight
with exhaustive modular determinants and random entry-flip cubes with exact
Bareiss determinants.  The optional tie output is independently verified and
stores the first nonzero frontier mask; reports retain the first 32 such masks
and explicitly flag truncation.

`fast_cube_batch.py` deterministically mixes top-single, pair-linkage,
neutral-plus-transverse, balanced-random, and destroy/repair supports.  It
deduplicates affine cubes, invokes only the audited evaluator, records
per-cube `/usr/bin/time` evidence, and atomically refreshes a resumable
aggregate report:

```sh
python3 research/fast_cube_batch.py \
  --output-dir runs/direct-search/fast-principal-cube/batch \
  --maximum-seconds 600

python3 research/fast_cube_tie_audit.py \
  --batch-dir runs/direct-search/fast-principal-cube/batch
```

`fast_cube_lnps.py` turns the per-weight records into a lineage-balanced
large-neighborhood beam.  It starts from H0/H1/H2 frontier representatives,
checks the two exact 12-flip bridges, constructs 27-entry supports from exact
one/pair ranks, and recenters with exact 6/9/12/15 overlaps.  Every admitted
beam matrix is independently Bareiss-verified; affine cubes and
sign-normalized states are separately deduplicated:

```sh
python3 research/fast_cube_lnps.py \
  --output-dir runs/direct-search/fast-principal-cube/lnps \
  --maximum-seconds 315
```

`fast_cube_partition32.py` exactly partitions one 32-entry affine cube into
32 disjoint 27-bit leaves.  Singular leaf bases are deterministically
rerooted within the same leaf and masks are XOR-mapped back to the global
32-bit assignment.  The report proves the `2^32` union, fingerprints the full
cube and every leaf, aggregates a globally exact top-K, and arena-verifies
all returned frontier ties. An explicit five-coordinate outer support permits
an audited full 32-entry endpoint connector:

```sh
python3 research/fast_cube_partition32.py \
  --output-dir runs/direct-search/fast-principal-cube/partition32

python3 research/fast_cube_partition32.py \
  --output-dir runs/direct-search/fast-principal-cube/connector32 \
  --start runs/start.matrix.txt \
  --endpoint runs/aligned-endpoint.matrix.txt \
  --calibrated-support runs/inner27.coords.txt \
  --outer-support runs/outer5.coords.txt \
  --bridge-size 32
```

The July 28 pilot evidence, tie-class audit, two completed 32-cubes, and
important assignment-visit/fingerprint caveats are summarized in
`research/FAST_CUBE_LNPS_20260728.md`.

`qd_cube_campaign.py` selects audited corrected-QD frontier outputs, restricts
supports to the dephased 22×22 core, and plans four exact 27-bit
adjugate/pair-informed supports per center plus optional 28–32-bit bridge
controls. Plans pin all runtime inputs, partition bridges into disjoint
27-bit tasks, and can be pooled safely by modulo shards:

```sh
python3 research/qd_cube_campaign.py plan \
  --output-dir runs/qd-selected-cubes

python3 research/qd_cube_campaign.py run \
  --output-dir runs/qd-selected-cubes \
  --shard-count 8 --shard-index 0 --resume

python3 research/qd_cube_campaign.py summarize \
  --output-dir runs/qd-selected-cubes
```

`affine_gf2_audit.py` maps full-matrix entry flips to exact dephased-core
GF(2) directions. It proves direction rank and detects affine intersections
or containment that equality-only cube fingerprints miss:

```sh
python3 research/affine_gf2_audit.py \
  --start runs/start.matrix.txt \
  --support runs/support.coords.txt \
  --label candidate-cube
```

The completed QD-selected tranche and a prior-corpus-disjoint 32-bit
connector are summarized in `research/QD_CUBE_CAMPAIGN_20260728.md`.

The recurrence follows the fast principal-minor method of Griffin and
Tsatsomeros, *Linear Algebra and its Applications* 419 (2006),
<https://doi.org/10.1016/j.laa.2006.04.008>.  The implementation here was
written for this repository; no archived third-party source was copied.

## Elite alignment and path relinking

`align_elite.py` heuristically minimizes Hamming distance between two factors
under signed row and column permutations. Alternating Hungarian assignments
make distant but equivalent-looking representations more useful as path
endpoints. It requires NumPy and SciPy, writes the aligned matrix and
transformation metadata atomically, and checks with Bareiss that the
determinant was preserved. New metadata binds the target, source, and aligned
output bytes by SHA-256:

```sh
python3 research/align_elite.py \
  --target runs/near.matrix.txt \
  --source runs/frontier.matrix.txt \
  --output runs/frontier-aligned.matrix.txt \
  --metadata runs/frontier-aligned.json \
  --seed 29640 --restarts 20000
```

A nonzero distance is not evidence of inequivalence; this is an optimization
helper, not a canonicalizer.

`path_relink.cpp` searches the interior of paths between repeated `--elite`
matrices. At each Hamming depth it ranks every one-step extension with inverse
updates, exactly checks the strongest `--exact-pool`, and retains an
exact-score beam plus a bounded diversity reserve. This reaches the middle of
large structured switches that a fixed-radius multiflip search cannot see:

```sh
c++ -std=c++20 -O3 -Wall -Wextra -Wshadow -Werror \
  research/path_relink.cpp -o build/research/path_relink

build/research/path_relink \
  --elite runs/near.matrix.txt \
  --elite runs/switch-1.matrix.txt \
  --elite runs/switch-2.matrix.txt \
  --output runs/path-best.matrix.txt \
  --log runs/path.jsonl --snapshot runs/path.snapshot.json \
  --seed 29630 --beam-width 1000 --exact-pool 4000 \
  --random-fraction 0.18 --seconds 3600
```

The parent-to-elite directions run before paths among non-parent elites.
`--max-pairs` can stop after that star; `--one-way` omits reverse paths.
Candidate promotion requires an exact Bareiss determinant. Matrix checkpoints,
the latest telemetry snapshot, and the append-only JSONL records are all
durably written. New telemetry records include the full elite paths, scores,
and collision-free row-major sign encodings, and output paths must be fresh.
Still run `./arena verify` before making a score claim.
This is a pruned heuristic search, not an exhaustive audit, and numerically
singular intermediate states are not expanded. Use immutable copies as elites:
an active search may replace its own output checkpoint while this campaign is
being prepared.

## Sharded radius-four screen

`radius4_screen.cpp` visits every four-entry perturbation assigned to its
shard. It derives the inverse from exact integer cofactors, evaluates the
rank-four determinant lemma in binary64, sends every prediction within
`--exact-margin` of the incumbent through Bareiss, and exactly checks the
retained `--top-pool` at the end:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Wshadow -Werror \
  research/radius4_screen.cpp -o build/research/radius4_screen

build/research/radius4_screen \
  --start runs/frontier.matrix.txt \
  --output runs/radius4-shard0.matrix.txt \
  --research-output runs/radius4-shard0-best-below.matrix.txt \
  --log runs/radius4-shard0.jsonl \
  --shard-count 4 --shard-index 0 \
  --top-pool 512 --exact-margin 1000000000000 \
  --seed 29660 --seconds 3600
```

The outermost flipped-entry index modulo `--shard-count` partitions the
`C(529,4) = 3,226,076,876` combinations without overlap. A shard is complete
only when its exact expected count matches `screened_combinations`; all shard
snapshots must be complete before claiming full numerical coverage.
New snapshots also embed the start path and a collision-free row-major sign
encoding. Older snapshots made before that field was introduced require a
separate digest-bound manifest before attributing their counts to a particular
start matrix.

This is deliberately labeled a **complete floating screen**, not an exact
radius-four audit. Calibration samples and the wide exact gate make numerical
misses observable and unlikely, but they are not a formal rounding-error
certificate. The global output never regresses; the separate research output
preserves the strongest exactly checked four-flip state in that shard.
Candidate, telemetry snapshot, and JSONL writes are durable, and any promoted
matrix still requires `./arena verify`. All output paths must be fresh; the
tool refuses to overwrite an earlier result or snapshot.

## Structured random seeds

Create a reproducible seed by replacing whole lines with fresh random sign
patterns:

```sh
python3 research/structured_seed.py \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/seed-23001.matrix.txt \
  --seed 23001 --lines 6 --orientation mixed
```

`rows` and `columns` replace exactly `--lines` lines of that orientation.
`mixed` splits the total between rows and columns, assigning the extra line to
rows, and therefore accepts totals from 2 through 46. The output is written
atomically and is suitable for `./arena verify` or as another search start.

To derive a seed from one of the 60 inequivalent order-24 Hadamard
representatives in the Mendeley archive, use:

```sh
python3 research/hadamard24_seed.py \
  --input /path/to/Hadamard_24.txt \
  --class 51 --delete-row 1 --delete-column 1 \
  --output runs/h24-class-51.matrix.txt
```

The parser first verifies `HH^T = 24I`; the emitted minor is written
atomically.

## Exact two-line audit

`pair_search.cpp` jointly optimizes two complete rows, then two complete
columns. For each pair it fixes one redundant sign and exhaustively evaluates
all `2^22` assignments; the other line is then chosen exactly from the
second-cofactor signs.

```sh
c++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  research/pair_search.cpp -o build/research/pair_search

build/research/pair_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/pair-audit.matrix.txt \
  --log runs/pair-audit.jsonl \
  --seed 2301 \
  --passes 1

./arena verify runs/pair-audit.matrix.txt
```

Use `--passes 1` for a completion-counted audit of all 506 unordered row and
column pairs. `--seconds N` is useful for exploratory repeated passes, but a
time-limited partial pass does not establish full neighborhood coverage.

For exploration beyond a line-local optimum, add `--kick-size K` with a time
budget. Each pass flips exactly `K` distinct entries of the global incumbent,
then accepts exact improvements in the two-line neighborhood while keeping the
best output checkpoint separate from the lower-scoring working state:

```sh
build/research/pair_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/pair-kick.matrix.txt \
  --research-output runs/pair-kick-best-below-frontier.matrix.txt \
  --log runs/pair-kick.jsonl \
  --seed 2401 \
  --kick-size 12 \
  --settle-passes 3 \
  --seconds 3600
```

`--output` is never replaced by a lower score. The optional, distinct
`--research-output` preserves the strongest working state still below the
global incumbent, which is useful for studying alternate basins. The tool
rejects aliases among its start, global output, research output, and log paths.
`--settle-passes N` lets an improving kicked basin complete multiple pair
sweeps before it is abandoned; a zero-move sweep still settles it early.

## Exact row-column cross audit

`cross_search.cpp` jointly replaces one complete row and one complete column.
For each of the 529 crosses it fixes the redundant global sign and exhausts
`2^21` assignments; the remaining signs are chosen exactly from the bilinear
cofactors.

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/cross_search.cpp -o build/research/cross_search

build/research/cross_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/cross-audit.matrix.txt \
  --log runs/cross-audit.jsonl \
  --seed 29101 --passes 1 --heartbeat-seconds 30
```

One completed pass is an exact audit of all
`529 * 2^21 = 1,109,393,408` assignments. With `--kick-size N`, the same
engine explores perturbed basins while retaining the best exact checkpoint.
Equal-score factors replace the incumbent only when they are not related by
row/column sign changes, keeping distinct cross neighborhoods in play without
cycling through sign-only variants.

## Exact Gram screens

`gram_edge_search.py` searches restricted defect graphs for a necessary Gram
condition before any decomposition attempt. It supports the published
reference toggle/swap neighborhoods, two-for-two reference swaps, and up to
four additions to the `K3 + 5 K4` Ehlich graph.

```sh
python3 research/gram_edge_search.py \
  --mode reference-double-swap \
  --shard-count 4 --shard-index 0 \
  --output runs/gram-double-swap-shard0.json
```

Determinants use an exact determinant lemma with independently Bareiss-checked
cache patterns. A surviving Gram would still require a separate
`G = AA^T`, `A in {-1,+1}^{23x23}` decomposition; the screen alone never
claims a matrix or an improved score.

`gram_radius3_orbits.cpp` exactly quotients the published graph's
delete-three/add-three neighborhood by its full 442,368-element automorphism
group. Build an immutable catalog first, then screen it unsharded or by
deterministic orbit-index shards:

```sh
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wshadow -Werror \
  research/gram_radius3_orbits.cpp \
  -o build/research/gram_radius3_orbits

build/research/gram_radius3_orbits --self-test

build/research/gram_radius3_orbits \
  --build-catalog --catalog runs/gram-radius3/catalog.tsv \
  --manifest runs/gram-radius3/catalog-manifest.json

build/research/gram_radius3_orbits \
  --screen --catalog runs/gram-radius3/catalog.tsv \
  --catalog-sha256 CATALOG_SHA256 \
  --output runs/gram-radius3/screen.json \
  --route-snapshot runs/gram-radius3/route.json
```

The catalog SHA is mandatory and screening reports distinguish assigned-shard
completion from whole-family completion. The retained complete result and
full artifact hashes are in `GRAM_RADIUS3_ORBITS_20260728.md`.

`gram_published_degree_slice.cpp` exactly enumerates the fixed ideal-base
slice obtained by adding 12 absent edges with added degrees `4,4,4` on the
published `K3` and `2^6,0^14` on the five-`K4` side. It quotients 3,488,400
labeled configurations to 20 exact base-automorphism orbits and audits overlap
with the connector and radius catalogs. None of its square determinants is
strictly above the frontier. See `GRAM_PUBLISHED_DEGREE_SLICE_20260728.md`.

`gram_radius4_basin.cpp` takes the strict-above-frontier representatives from
the frozen radius-3 orbit catalog and exactly screens every true outward
radius-4 child.
Its primary rank-4 parent update and independent rank-16 reference update must
agree; the catalog SHA is pinned in the binary. The retained run screened
`16,858,380` parent-child transitions, routed four `Aut(B0)` square/PD
orbits, and rejected all four with exact shell certificates. This is a
frontier-parent basin closure, not a full radius-4 enumeration. See
`GRAM_RADIUS4_BASIN_20260728.md`.

## Stochastic Gram graph search

`gram_tabu.cpp` searches normalized graphs `G = 24I - J + 4A` in hill, tabu,
or annealing mode. Floating-point Cholesky and determinant-lemma values choose
moves; every one-edge neighbor is screened by exact CRT reconstruction.
Stored hits must have an exact square determinant, a square root divisible by
`2^22`, and exact positive leading minors.

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/gram_tabu.cpp -o build/research/gram_tabu

build/research/gram_tabu \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/gram-tabu.json \
  --mode tabu --seed 29303 --seconds 3600 \
  --kick-size 11 --tabu-tenure 17 \
  --max-stored-hits 20000
```

The checkpoint is atomic, bounded by `--max-stored-hits`, and records both
stored-hit and omitted-observation counters. `SIGINT`/`SIGTERM`, wall limits,
and iteration limits write a final `complete: true` snapshot. A hit is only a
Gram candidate: apply the exact obstruction filters below, then require an
exact sign decomposition before making any matrix claim.

## Multilevel Gram search

`gram_multilevel.cpp` covers the nearby Gram space missed by the binary graph
model: diagonal entries are `23`, while off-diagonal entries may be
`-5`, `-1`, or `3`, with optional `7` entries. It row-switches a supplied
sign matrix into that normalization, keeps every accepted search state exactly
positive definite, and exact-screens every adjacent-level neighbor. Only
determinant squares strictly above the frontier are stored as hits.

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/gram_multilevel.cpp -o build/research/gram_multilevel

build/research/gram_multilevel \
  --start runs/near-frontier.matrix.txt \
  --output runs/gram-multilevel.json \
  --mode tabu --seed 29611 --seconds 3600 \
  --min-minus5 1 --max-minus5 4 \
  --kick-size 11 --tabu-tenure 19
```

The `-5` bounds matter: unconstrained determinant ascent normally removes the
rare `-5` edge and falls back into the binary search space. Snapshots bind the
raw start matrix by SHA-256 and atomically record parameters, exact counters,
stored hits, capacity omissions, and termination. Both exact obstruction
filters below accept the binary and multilevel snapshot schemas. The current
framed decomposer accepts only binary snapshots, so a passing multilevel hit
still needs a multilevel-capable exact decomposition step.
The engine has no resume mode and therefore refuses to replace an existing
snapshot at startup.
With a positive `-5` minimum, `--relocate-minus5` adds an atomic two-edge move
that relocates one `-5` edge while preserving the current count. This works
both on fixed-count strata and inside a bounded count range. The corresponding
`--min-plus7`, `--max-plus7`, and `--relocate-plus7` options provide the same
bounded relocation for `7` defects. Four-level runs use a larger exact CRT
modulus whose product exceeds twice the `1607^(23/2)` row-norm determinant
bound.

## Exact rational obstruction

Before spending decomposition time, `gram_hasse.py` can exactly reject a Gram
candidate that is not rationally congruent to the identity:

```sh
python3 research/gram_hasse.py \
  --snapshot runs/direct-search/gram-tabu-long-29302.json \
  --hit-index 0 --output runs/gram-hit-0-hasse.json
```

The tool reconstructs any supported normalized Gram schema, checks its
determinant and positive definiteness, factors the exact square root with a
bounded deterministic Pollard--Brent implementation, and evaluates the Hasse
invariant at `2` and every prime dividing the determinant. `bad_primes` is an
exact rational obstruction to `G=RR^T`. An empty list is only “no
obstruction”: it neither constructs a factor nor implies a `{-1,+1}`
decomposition. Repeat
`--hit-index` to select several stored hits, or use `--all-hits` to scan every
hit retained in that snapshot. Capacity-omitted and deduplicated observations
are not silently classified; the report binds the snapshot SHA-256 and copies
its completion, capacity, and omission counters. Apply the separate exact
column-shell finite-field gate to survivors before spending time on framed
decomposition.

Factoring is deliberately limited to 64-bit roots (enough for order 23) and a
fixed work budget. Budget exhaustion produces `status: error`, never a false
rejection.

`gram_shell_filter.cpp` applies the stronger sign-specific gate. If
`G=RR^T`, every column sign vector `s` of `R` satisfies
`s^T G^-1 s = 1`, and those column outer products sum to `G`. The tool
computes `G^-1=M/d` exactly from CRT cofactors, exhausts all `2^22` sign
vectors up to global negation, and tests over `F_3` whether `G` belongs to
the span of their outer products:

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/gram_shell_filter.cpp -o build/research/gram_shell_filter

build/research/gram_shell_filter \
  --snapshot runs/direct-search/gram-tabu-long-29302.json \
  --hit-index 166 --hit-index 169 \
  --output runs/gram-shell-filter.json
```

`--hit-index` is repeatable and `--all-hits` filters all stored hits in a
snapshot.
Snapshot mode normally requires a square root strictly above the frontier.
`--allow-subfrontier-research` explicitly relaxes only that root threshold
for tied or lower research Grams; all square, divisibility, positive-definite,
snapshot-integrity, and shell checks remain exact, and the output records that
the opt-in mode was used.
`--matrix references/orrick-et-al-2003/matrix.txt` exercises the gate on a
known factor; matrix mode first row-switches any compatible Gram into the
supported `{-5,-1,3,7}` normalization. A rejection records the exact shell
size, span and augmented ranks, and a sparse mod-3 dual certificate that
annihilates every shell outer product but not `G`. Snapshot determinant, square,
divisibility, edge uniqueness, normalization, and positive definiteness are
independently recomputed, and `snapshot_sha256` binds the certificate to the
exact input bytes. Matrix-mode results likewise record the raw matrix SHA-256,
and each result names its own normalization. A pass remains only a
necessary-condition pass.

For a shell pass, `--include-shell-vectors` adds every normalized shell sign
vector to the JSON as a compact 23-bit mask. `gram_shell_milp.py` then searches
for nonnegative integer column multiplicities whose outer products sum exactly
to `G`:

```sh
build/research/gram_shell_filter \
  --matrix references/orrick-et-al-2003/matrix.txt \
  --include-shell-vectors --output runs/reference-shell-columns.json

python3 research/gram_shell_milp.py \
  --shell-report runs/reference-shell-columns.json \
  --output runs/reference-shell-factor.matrix.txt \
  --metadata runs/reference-shell-factor.json \
  --time-limit 300
```

The MILP solver is a discovery mechanism, not a proof of infeasibility. A
candidate factor is retained only after exact outer-product reconstruction,
exact Gram equality, and Bareiss determinant checks. `--seed N` uses a
reproducible random linear objective to sample alternate factors of the same
Gram. For a nonsingular target, `--binary --positive-objective
--mip-rel-gap 1` is the faster feasible-first formulation: normalized columns
cannot repeat, and the nonzero gap is used only to obtain a feasible proposal,
never to claim objective optimality. Repeatable `--exclude-factor` options add
exact raw-subset no-good cuts.

`gram_shell_orbit_cpsat.py` specializes the same Boolean factor problem to the
published frontier Gram. It derives the exact `Aut(G)` shell orbits, proves the
forced `3,6,6,8` orbit incidences, and fixes one of four canonical triples in
the size-six orbit before calling CP-SAT. Solver `UNKNOWN` remains only a
bounded search result. The retained campaign found a new local H/HT factor
class and is documented in `FRONTIER_FACTOR_CLASS_EXPANSION_20260728.md`.
The subsequent two-arm search and exhaustive `2^32` connector from that class
are documented in `EB138A_SEARCH_20260728.md`.

The later novelty-aware factor harvest and two sparse integer identities close
the other three size-six triple orbits exactly for this fixed Gram, replacing
their earlier `UNKNOWN` statuses with checkable linear contradictions. The
pilot, proof boundary, standalone checker, and local-versus-literature claims
are documented in
[`FRONTIER_PORTAL_HARVEST_20260729.md`](FRONTIER_PORTAL_HARVEST_20260729.md).

The full 2026-07-29 portal campaign synthesis, including the two new local
H/HT classes, exact alignment and radius-three closures, pairwise portal
geometry, four exhaustive 32-coordinate connector cubes, and completed
core/fiber follow-ups, is documented with terminal statistics and provenance
hashes in
[`FRONTIER_PORTAL_CAMPAIGN_20260729.md`](FRONTIER_PORTAL_CAMPAIGN_20260729.md).

Always pass an emitted matrix through `./arena verify` before promotion.

## Randomized Gram decomposition

`gram_decompose.cpp` implements the framed row-decomposition algorithm of
Brent, Orrick, Osborn, and Zimmermann. It derives `G = AA^T` from a supplied
factor, quotients signed column permutations with prefix frames, and solves the
bounded integer constraints for each new row exactly. Random restarts retain
one child and, by default, a second child with probability `0.3`.

Because this arena already has a factor, the default search keeps its framed
path as an anchor and samples sibling subtrees. This guarantees a known leaf
without weakening the exact checks on generated factors. Use `--no-anchor` for
the literature's pure randomized tree search, or add more known factors with
repeatable `--anchor-factor` options. Supplied factors seed deduplication and
do not consume `--max-solutions`; only newly generated column-canonical
factors are emitted.

```sh
c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/gram_decompose.cpp -o build/research/gram_decompose

build/research/gram_decompose \
  --input references/orrick-et-al-2003/matrix.txt \
  --output-dir runs/gram-factors \
  --log runs/gram-factors.jsonl \
  --seed 23003 --seconds 60 --max-solutions 10
```

To attack a qualified, previously unfactored survivor from `gram_tabu.cpp`,
select its zero-based hit directly from the exact checkpoint:

```sh
build/research/gram_decompose \
  --gram-snapshot runs/direct-search/gram-tabu-smoke.json \
  --hit-index 0 \
  --output-dir runs/gram-hit-0-factors \
  --log runs/gram-hit-0-decompose.jsonl \
  --seed 23008 --seconds 60 --max-solutions 10
```

Snapshot mode requires the exact `G=24I-J+4A` normalization and a hit marked
qualified. It reconstructs the 23 by 23 Gram matrix from the one-based edge
list, then independently checks edge uniqueness, the claimed square,
`2^22` divisibility, the exact CRT determinant, and every positive leading
principal minor. It never trusts the snapshot's numerical claims alone.

Every emitted matrix is independently checked with integer arithmetic for
`RR^T = G` and the target Bareiss determinant, then written atomically under
its SHA-256 column-canonical hash. JSONL records include restarts, tree depth,
bounded assignments, feasible/sample counts, exact checks, and heartbeats.
This hash removes signed column permutations only; it does not establish full
Hadamard inequivalence between emitted factors. The tool uses POSIX atomic
writes and signed 128-bit exact arithmetic and is not intended as a portable
MSVC utility; concurrent runs should use distinct output directories and log
files.
The algorithm is described in Section 4.1 of
https://arxiv.org/abs/1112.4160 and randomized in Section 2 of
https://arxiv.org/abs/1112.4671.

## Hadamard descendant campaigns

The 60 order-24 representatives yield 1,106 non-isomorphic Hadamard
`2-(23,11,5)` designs after dephasing, so the default `(1,1)` minor from each
class is only a small slice of the extendable order-23 family.
`h24_deletion_campaign.py` checks the Mendeley representatives against the
local Spence design catalog, enumerates all 576 row/column pivots for selected
classes, runs a short reactive search from each minor, independently
Bareiss-checks every retained score, and writes bounded atomic checkpoints.

```sh
python3 research/h24_deletion_campaign.py \
  --classes 14,42,51 \
  --seconds-per-minor 0.25 \
  --seed-base 32000 \
  --top-count 30 \
  --output-dir runs/h24-deletion-14-42-51
```

Use sharding for a full 60-class campaign. Pivots are assigned a seed from
their class, row, and column independently of the shard. See
`research/H24_STRUCTURED_20260728.md` for the coverage audit and exact
artifacts.

## Exact order-22 bordering

Current audited plateau and fixed-border results are summarized in
[`ORDER22_PLATEAU_ATLAS_20260729.md`](ORDER22_PLATEAU_ATLAS_20260729.md).
The seeded order-23 follow-up waves and exact local postpasses are summarized
in [`ORDER22_SEEDED_ORDER23_SEARCH_20260729.md`](ORDER22_SEEDED_ORDER23_SEARCH_20260729.md).

`corpus_matrix.py` extracts a strict sign matrix from the local Mendeley
maximal-determinant corpus. `order22_border.cpp` then fixes that 22 by 22 core
and exhausts all `2^21` border columns. For each column `x`, it chooses the
globally optimal border row directly from `adj(B)x`; every score is an exact
integer and every promoted 23 by 23 matrix is Bareiss-checked.

```sh
python3 research/corpus_matrix.py \
  --order 22 --output runs/order22.matrix.txt

c++ -std=c++20 -O3 -march=native -Wall -Wextra -Werror -pedantic \
  research/order22_border.cpp -o build/research/order22_border

build/research/order22_border \
  --start runs/order22.matrix.txt \
  --output runs/order22-border.matrix.txt \
  --log runs/order22-border.jsonl
```

This is exhaustive only for the specified fixed core; it is not an order-23
global statement.

## Signed-perfect-matching cycle pilot

Dephased order-23 matrices can be read as bipartite graphs whose determinant
is the signed perfect-matching imbalance. The exact bounded `k=3,5`
same-sign alternating-cycle pilot, including its negative result and receipt
audit, is summarized in
[`MATCHING_CYCLE_PILOT_20260729.md`](MATCHING_CYCLE_PILOT_20260729.md).

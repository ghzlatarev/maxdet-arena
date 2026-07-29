# Exact Hamming-code clique LNS v2 (2026-07-29)

## Hypothesis

After signed row normalization, every off-diagonal order-23 Gram entry is
`3 mod 4`. Anchoring the first row at all `+1` turns the high-quality Gram
alphabet

```text
inner product    3   -1   -5
Hamming distance 10   12   14
```

into a binary-code constraint. This suggests replacing a large set of rows at
once as an exact clique-completion problem rather than perturbing matrix
entries independently.

`research/hamming_clique_lns.cpp` implements that representation as a strict,
time-bounded exact repair engine. In each iteration it:

1. chooses `t` non-anchor rows to destroy;
2. exhaustively tests all `2^23` oriented row words against the fixed rows;
3. builds the exact bitset compatibility graph on the survivors;
4. uses a proper greedy coloring to reject candidate subgraphs that cannot
   contain the required clique;
5. enumerates the remaining `t`-cliques under an exact fixed-Gram Schur
   residual bound; and
6. ranks every unpruned completion with a `cpp_int` Bareiss determinant.

The original v1 engine used the sound but loose bound

```text
det(R R^T) * 23^(23-|R|) <= incumbent^2.
```

This is the basic Hadamard--Fischer bound: each remaining Schur diagonal is at
most the original row norm squared, `23`. Version 2 replaces it with the
strictly more informative exact bound below.

## Exact v2 Schur bound

Let `F` contain the `23-t` fixed rows, and define

```text
A = F F^T
D = det(A)
J = adj(A).
```

For candidate row `x_i`, let `b_i = F x_i^T`. The dense integer Schur kernel
is

```text
N_ij = D <x_i,x_j> - b_i^T J b_j.
```

For a chosen candidate set `Q` of size `k`, block elimination gives the exact
identity

```text
det(G_Q) = det(N_Q) / D^(k-1),
```

where `G_Q` is the Gram matrix of the fixed rows followed by `Q`. To reduce
every later operand, v2 computes

```text
g = gcd(D, every N_ij)
N = g K
D = g d.
```

It follows that

```text
det(G_Q) = g det(K_Q) / d^(k-1).
```

For `k=0`, read this directly as `det(G_empty)=D`; the integral pruning
formula below handles the empty choice without using a negative power.

Now let `r=t-k`, let `C` be the current DFS candidate tail, and let `P_r` be
the product of the `r` largest diagonals `K_ii` for vertices in `C`.
The fixed-row residual of candidate `i` is `K_ii/d`. Conditioning on the rows
already in `Q` can only decrease that residual. Hadamard--Fischer therefore
gives

```text
det(full Gram) <= g det(K_Q) P_r / d^(t-1).
```

For incumbent absolute determinant `T`, v2 prunes exactly when

```text
g det(K_Q) P_r <= T^2 d^(t-1).
```

All factors and comparisons are `cpp_int`; floating-point arithmetic is not
used to prune. For the empty choice, `det(K_empty)=1`. For a full choice,
`P_0=1`, so the endpoint comparison is exact. A zero principal determinant
means the selected rows are dependent and every extension is singular.

The small principal determinant is updated by exact bordered identities. Its
adjugate is materialized lazily only when the DFS actually descends. Before
that determinant work, a proper coloring partitions the current candidate
subgraph into independent color classes. If it uses fewer than `r` colors,
the subgraph cannot contain an `r`-clique and is rejected soundly.

## Build and run

On the current macOS development host, Boost headers are under Homebrew:

```sh
clang++ -std=c++20 -O3 -I/opt/homebrew/include \
  -Wall -Wextra -Wshadow -Wpedantic -Werror \
  research/hamming_clique_lns.cpp \
  -o build/research/hamming_clique_lns

build/research/hamming_clique_lns \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/hamming-clique/best.matrix.txt \
  --log runs/hamming-clique/run.jsonl \
  --summary runs/hamming-clique/summary.json \
  --seed 2306 --seconds 300 \
  --destroy-min 4 --destroy-max 7 \
  --destroy-without-replacement
```

Use `--allow-distance-14` only for the separate arm that admits normalized
Gram entry `-5`. The default arm is the tighter `{10,12}` code containing the
published comparison matrix.

Use `--transpose-start` to run the identical row-repair search on the input's
column neighborhoods without materializing a second matrix file.

`--destroy-without-replacement` ranks the masks of each requested destroy
size over the 22 non-anchor rows and consumes each mask at most once for the
current parent. `--destroy-shard-count N --destroy-shard-index I` assigns mask
rank `r` to shard `r mod N`, so the shards are disjoint and their union is the
complete mask family for a shared parent. Within a shard, a pinned SplitMix64
key derived from the seed, parent generation, destroy size, and mask gives a
deterministic order. The scheduler samples uniformly among destroy sizes that
still have masks, rather than weighting sizes by their numbers of masks.

A strict improvement increments the parent generation and deterministically
resets the schedules because the same mask is then a different neighborhood.
Mask, rank, parent generation, and shard coordinates are recorded per
iteration. Exhausting all selected masks ends with
`destroy_schedule_exhausted`. Supplying either shard option implicitly enables
without-replacement scheduling.

The output, JSONL log, and JSON summary must be new, mutually distinct paths.
The executable refuses to overwrite any of them. It installs and updates each
owned artifact by same-directory atomic replacement and polls both the wall
deadline and `SIGINT`/`SIGTERM`.

## Validation

The retained validation commands were:

```sh
build/research/hamming_clique_lns --self-test

clang++ -std=c++20 -O1 -g -I/opt/homebrew/include \
  -Wall -Wextra -Wshadow -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  research/hamming_clique_lns.cpp \
  -o /tmp/hamming_clique_lns_asan
/tmp/hamming_clique_lns_asan --self-test
```

The self-test independently checks:

- switch normalization and the known `det(J-2I)` value;
- complete compatible-word enumeration on small cases;
- every bitset edge, proper color class, and all size-three cliques against
  brute combinations;
- exact matrix/Gram determinant identity;
- the original Hadamard--Fischer bound on deterministic random small
  matrices;
- the Schur/gcd determinant identity, lazy bordered determinants and
  adjugates against independent Bareiss/cofactor calculations, and the
  top-`r` residual bound against brute completions, including empty, full, and
  singular selections; and
- destroy-mask rank/unrank, no repetition, shard disjointness and complete
  union, seed-independent membership, deterministic order and parent reset,
  schedule exhaustion, and invalid shard arguments.

An earlier v1 ASan/UBSan one-iteration search smoke also completed cleanly. A
v1 release smoke destroying six rows of the published matrix examined all
`8,388,608` oriented words, retained 59 candidates, checked 1,711 graph
pairs, traversed 8,579 clique nodes, and completed in 0.46 seconds. Its
emitted matrix independently passed `./arena verify` at
`2,779,447,296,000,000`. This was a correctness/performance smoke, not a
search campaign, and it found no strict improvement.

The distance-14 arm was separately smoked on a normalized seed containing one
`-5` Gram pair. The default arm correctly rejected that seed; the explicit
distance-14 arm completed and its emitted matrix independently passed the
arena verifier.

An independent proof-and-implementation audit of source SHA-256
`97eb017e4f9df338b6f1594b862cf48e70140531030aaca56e499add24f4076d`
checked 32,314 randomized or exhaustive cases. Its scope included the
Schur/gcd identity, the actual current DFS candidate tail used by the top-`r`
bound, empty and full choices, singular branches, signed gcd inputs,
incumbent updates, proper coloring, and complete/disjoint sharding per parent
generation. Strict compilation, the built-in self-test, ASan/UBSan, and a real
frontier `t=8` repair also passed. The audited v2 output was byte-identical to
the matched v1 output. This is an implementation soundness audit, not a claim
about unsearched neighborhoods or global optimality.

## Matched v2 benchmark

A matched one-iteration benchmark used the published frontier, seed `42301`,
destroy size exactly eight, the default `{10,12}` alphabet, and the same
350-vertex pool. The original and lazy-bordered columns are v1 controls; the
v2 column is the finalized Schur, coloring, and scheduling-capable engine.

| Measure | Original v1 | Lazy-bordered v1 | Final v2 |
| --- | ---: | ---: | ---: |
| Elapsed seconds | 81.2007 | 27.2763 | 0.079474 |
| Exhaustive 23-bit word tests | 8,388,608 | 8,388,608 | 8,388,608 |
| Candidate vertices | 350 | 350 | 350 |
| Clique branch nodes | 3,446,739 | 3,446,739 | 27,262 |
| Cardinality prunes | 2,656,298 | 2,656,298 | 5,330 |
| Proper-color checks | — | — | 22,245 |
| Proper-color prunes | — | — | 17,326 |
| Exact determinant-bound checks | 1,160,331 | 1,160,331 | 4,921 |
| Determinant-bound prunes | 790,439 | 790,439 | 4,605 |
| Full leaves scored | 0 | 0 | 0 |
| Improvements | 0 | 0 | 0 |

All three runs completed the same exact destroyed-row neighborhood, retained
the frontier determinant, and emitted byte-identical matrices. The timings
are single-host engineering measurements, not a general performance
guarantee. No optimized multi-arm v2 search result is inferred from this
benchmark.

## Retained depth calibration

`runs/hamming-clique-schur-depth-calibration-20260729` retained one fixed-seed
frontier arm at each of destroy sizes 9, 10, and 11:

| Destroy size | Pool | Clique nodes | Elapsed seconds | Result |
| ---: | ---: | ---: | ---: | --- |
| 9 | 774 | 144,146 | 0.613866 | 1 exact neighborhood completed |
| 10 | 1,430 | 7,372,965 | 25.324854 | 1 exact neighborhood completed |
| 11 | 2,715 | 0 | 0.071033 | checked dense-kernel safety abort |

The size-11 pool required more than the allowed 4,000,000 dense Schur-kernel
entries, so the engine aborted explicitly without truncating the pool or
searching the clique neighborhood. All three emitted matrices independently
passed `./arena verify` at the frontier determinant. Only the completed
size-9 and size-10 runs are exact neighborhood closures, each for its one
recorded destroyed-row set and the `{10,12}` alphabet. These three isolated
measurements do not establish a general timing curve.

## First optimized v2 campaign

`runs/hamming-clique-schur-portal-wave-20260729` ran eight 600-second arms at
destroy size exactly eight. Four disjoint-shard arms covered row and
transposed neighborhoods of the published frontier; four unsharded arms used
the two locally distinct frontier-equivalent H/HT factors `de764` and `1e4b`
in both orientations.

| Measure | Combined |
| --- | ---: |
| Completed exact neighborhoods | 53,800 |
| Deadline-interrupted neighborhoods | 7 |
| Candidate words examined | 451,340,779,524 |
| Completed full candidate pools | 53,803 |
| Maximum pool size | 710 |
| Clique branch nodes | 551,659,248 |
| Proper-color prunes | 295,612,308 |
| Exact Schur-bound prunes | 151,999,647 |
| Full leaves requiring matrix scoring | 0 |
| Improvements | 0 |

All eight final matrices independently passed `./arena verify` at
`2,779,447,296,000,000`. The seven deadline-interrupted iterations are
excluded from the 53,800 closures; three of them completed their full pool
before stopping later in the graph search. The 451,340,779,524 count includes
all exact word tests executed, including partial interrupted pools.

Each completed iteration closes only its recorded eight-row or eight-column
set for the `{10,12}` alphabet. The frontier shards stayed at their common
unchanged parent and therefore remained disjoint, but they covered only the
recorded shard prefixes before the deadline. This is not a closure of every
size-eight subset or a statement about global optimality.

## Optimized size-nine campaign

`runs/hamming-clique-schur-t9-portal-wave-20260729` applied the same eight-arm
frontier/portal and row/transposed layout for 900 seconds per arm at destroy
size exactly nine.

| Measure | Combined |
| --- | ---: |
| Completed exact neighborhoods | 16,902 |
| Deadline-interrupted neighborhoods | 8 |
| Candidate words examined | 141,842,649,090 |
| Completed full candidate pools | 16,908 |
| Maximum pool size | 1,183 |
| Clique branch nodes | 2,340,939,470 |
| Proper-color prunes | 1,227,908,000 |
| Exact Schur-bound prunes | 741,823,879 |
| Full leaves requiring matrix scoring | 0 |
| Improvements | 0 |

All eight final matrices independently passed `./arena verify` at the
frontier. The eight interrupted iterations are excluded from the 16,902
closures; six completed their full pool before stopping later. The word-test
count includes partial interrupted pools.

Across the optimized size-eight and size-nine waves, 70,702 exact recorded
neighborhoods were completed, 593,183,428,614 candidate words were examined,
and 2,892,598,718 clique nodes were traversed. These are additive executed-work
counts, not a claim that the underlying global search regions are disjoint.
The claim boundary remains each recorded destroyed set, orientation, parent
generation, and `{10,12}` alphabet.

## Sparse size-eleven boundary crossing

Source SHA-256
`ebe4731fe25bf26e6499824906fb7b066ae926f5fb980d7cc941f3f788149f61`
stores exact Schur diagonals plus compatibility-edge entries rather than a
dense candidate-square matrix. Every principal kernel reached by clique DFS
contains only those entries, so taking `gcd(D, diagonals, edge entries)`
preserves the same determinant identity. Off-edge access is rejected, and the
allocation guard is checked in bytes without truncating candidates.

On the retained size-eleven fixture, the sparse kernel held 2,715 candidates
and 1,795,185 compatibility edges. Exact construction took 0.420 seconds; the
planned kernel/workspace was 57.34 MiB versus 224.95 MiB for the old dense
integer payload. The retained three-second full-run smoke crossed the former
refusal boundary at 80.73 MiB peak RSS, traversed 591,872 nodes, made 294,298
bound checks, pruned 291,191 of them, and found no leaf or improvement before
its deadline. The output remained at the frontier.

This smoke validates the implementation boundary and memory reduction; it is
not a completed size-eleven neighborhood and therefore adds no closure.
Strict compilation, self-tests, ASan/UBSan, and a matched size-eight replay
passed. The replay was byte-identical with all counters unchanged. Two
independent source audits found no soundness defect; their comparisons also
covered additional size-eight, size-nine, and distance-14 neighborhoods.

The checked 128 MiB guard budgets the explicitly modeled kernel arrays and
workspace, not total process RSS. Allocator overhead, multiprecision limbs,
and the already-built graph sit outside that estimate. Kernel construction
also does not poll the wall deadline. These are resource/latency caveats, not
approximations or changes to the exact pruning proof.

## First exact campaigns

Four five-minute row-repair arms started from the frontier, c1, c4, and the
high-`det(M)` core-quotient-632,590,000 basin. A matched four-arm campaign used
`--transpose-start` to repair column neighborhoods.

These campaign results were produced by the pre-v2 engine and are retained
unchanged as the historical search record.

| Measure | Row wave | Transposed wave | Combined |
| --- | ---: | ---: | ---: |
| Completed exact neighborhoods | 90 | 69 | 159 |
| Deadline-interrupted neighborhoods | 4 | 4 | 8 |
| Exhaustive 23-bit word tests | 788,529,152 | 612,368,384 | 1,400,897,536 |
| Clique branch nodes | 39,308,479 | 40,292,693 | 79,601,172 |
| Exact Fischer bound checks | 14,634,084 | 14,314,966 | 28,949,050 |
| Bound-pruned branches | 9,506,429 | 8,926,462 | 18,432,891 |
| Full leaves requiring matrix scoring | 0 | 0 | 0 |
| Improvements | 0 | 0 | 0 |

Every completed destroyed-line set was closed exactly for distance alphabet
`{10,12}`. Eight in-progress neighborhoods stopped at their wall deadlines
and are not counted as closures. All eight final arm outputs independently
passed `./arena verify`.

The result is useful negative local evidence, not a radius-eight closure:
destroyed subsets were sampled stochastically, and one difficult eight-line
repair can consume much of an arm. Future campaigns should stratify destroy
size instead of drawing uniformly from 4 through 8; the v2 no-replacement
scheduler now supports deterministic stratification and sharding.

## Limitations and claim boundary

- The candidate pool is exact or the iteration aborts; it is never truncated.
  `--pool-limit` is therefore a memory-safety refusal threshold, not a search
  approximation.
- The current edge-only exact Schur kernel stores every diagonal and
  compatibility edge. A checked byte guard rejects allocations that exceed
  its configured safety boundary; this remains an explicit abort, not
  candidate truncation or approximate search.
- A deadline can stop pool generation, graph construction, or clique repair.
  The summary distinguishes started, completed, and interrupted iterations.
- A completed iteration closes only its recorded destroyed row set in the
  selected distance alphabet. It says nothing about other row subsets,
  alphabets, switching classes, or global optimality.
- Without-replacement coverage is relative to one parent. If an improvement
  changes that parent, the generation and schedule reset; independently
  improving shards no longer describe the same family of neighborhoods.
- The graph search uses cardinality, proper-color feasibility, and exact
  Schur-residual Hadamard--Fischer pruning. Performance outside the measured
  frontier `t=8` pool, especially for much larger pools or the distance-14
  alphabet, remains empirical.
- Each arm replaces one orientation's rows only. Use `--transpose-start` in a
  separate arm to expose the corresponding column neighborhoods.
- The optimized campaigns still cover finite recorded prefixes of destroy
  schedules; they do not close a full destroy radius.

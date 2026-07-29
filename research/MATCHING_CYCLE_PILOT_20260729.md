# Signed-perfect-matching cycle pilot

## Representation and theorem

Dephase an order-23 sign matrix `A` so its first row and column are `+1`,
delete them, and set

```text
B[i,j] = (1 - A[i+1,j+1]) / 2.
```

Subtracting the first row of `A` from every other row gives the exact identity

```text
det(A) = (-2)^22 det(B) = 2^22 det(B).
```

Regard the binary matrix `B` as the adjacency matrix of a bipartite graph.
The determinant expansion says

```text
det(B) = (# even perfect matchings) - (# odd perfect matchings).
```

Fix a perfect matching `M` with sign `sign(det(B))` and permute its columns
onto the diagonal. If another matching differs from `M` on one alternating
cycle containing `k` matched pairs, its relative permutation is a `k`-cycle,
so its sign relative to `M` is

```text
(-1)^(k - 1).
```

Thus an odd alternating cycle gives a same-sign matching. This suggests two
concrete moves when all off-diagonal cycle edges are absent:

- **completion:** add the odd-cycle edges, retaining `M` and creating the
  same-sign alternative;
- **exchange:** delete the affected edges of `M` and add the odd-cycle edges,
  preserving every row/column degree while replacing `M` by a same-sign
  matching.

This theorem controls the sign of that one matching contribution. It does not
control the background matchings created or destroyed by the move, so an exact
experiment is required.

## Non-duplication boundary

The repository already implements the immediate maximum-volume-simplex
consequence: optimizing one cube vertex is the existing cofactor-sign
whole-line replacement. The research inventory contained no dedicated
perfect-matching-parity cycle generator. Generic affine-cube tools could in
principle evaluate one of these supports, so this is a new directed pilot in
the local method inventory, not a mathematical novelty claim.

## Bounded pilot

[`matching_cycle_search.py`](matching_cycle_search.py) used the fixed
first-row/first-column dephasing and one deterministic majority-sign matching
for each of:

- the published frontier matrix;
- all four retained order-23 outputs from the component5 campaign.

It enumerated directed absent cycles modulo cyclic rotation for `k=3,5` and
tested both moves. There were 110,131 cycle supports and 220,262 exact
candidates:

| input | source `abs(det)` | absent cycles | best retained `abs(det)` |
| --- | ---: | ---: | ---: |
| frontier | 2,779,447,296,000,000 | 20,364 | 2,458,818,969,600,000 |
| component5 arm0 | 2,687,827,968,000,000 | 23,425 | 2,346,172,022,784,000 |
| component5 arm1 | 2,687,827,968,000,000 | 21,830 | 2,419,200,098,304,000 |
| component5 arm2 | 2,687,827,968,000,000 | 19,986 | 2,408,448,000,000,000 |
| component5 arm3 | 2,621,033,676,800,000 | 24,526 | 2,341,416,730,624,000 |

No enumerated candidate improved its own source, and none beat the frontier.
The strongest frontier-rooted candidate was still
`320,628,326,400,000` below the frontier.

All candidate scores were evaluated with the exact adjugate form of the
determinant lemma. The best matrix in each of the 20
input/length/move families was then independently recomputed with Bareiss on
both its `22 x 22` binary core and reconstructed `23 x 23` sign matrix. All 20
retained matrices subsequently passed `./arena verify`, and the finalized
report binds every receipt and repeats verification independently.

The concise frozen manifest is
[`MATCHING_CYCLE_PILOT_20260729.json`](MATCHING_CYCLE_PILOT_20260729.json).
The full local report is
`runs/direct-search/matching-cycle-pilot-20260729/report.json`.

## Reproduce

```sh
python3 research/matching_cycle_search.py --self-test

python3 research/matching_cycle_search.py \
  --input frontier=references/orrick-et-al-2003/matrix.txt \
  --input component5-arm0=runs/direct-search/order22-component5-order23-search-20260729/arm0/best.matrix.txt \
  --input component5-arm1=runs/direct-search/order22-component5-order23-search-20260729/arm1/best.matrix.txt \
  --input component5-arm2=runs/direct-search/order22-component5-order23-search-20260729/arm2/best.matrix.txt \
  --input component5-arm3=runs/direct-search/order22-component5-order23-search-20260729/arm3/best.matrix.txt \
  --cycle-lengths 3,5 \
  --output-dir runs/direct-search/matching-cycle-pilot-20260729

for matrix_path in runs/direct-search/matching-cycle-pilot-20260729/*.matrix.txt
do
  ./arena verify "$matrix_path" \
    --json "${matrix_path%.matrix.txt}.receipt.json" --quiet
done

python3 research/matching_cycle_search.py \
  --finalize-report \
  runs/direct-search/matching-cycle-pilot-20260729/report.json
```

## Claim boundary

The negative result covers only absent `k=3,5` cycles at one fixed dephasing
and one deterministic majority-sign matching per input. It is not an audit of
all pivots, all perfect matchings, all alternating cycles, arbitrary cycle
edge subsets, or post-cycle local ascent. It proves neither local nor global
optimality. The useful negative evidence is narrower: naive single-matching
odd-cycle completion and degree-preserving exchange were uniformly poor on
these five audited high-determinant centers.

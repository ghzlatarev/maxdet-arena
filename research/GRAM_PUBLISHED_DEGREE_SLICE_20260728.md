# Exact published-degree Gram slice

**Run date:** 2026-07-28

## Scope

This audit fixes the ideal defect graph

```text
K3 disjoint-union 5 K4
```

and adds exactly 12 previously absent edges from the `K3` to the 20 vertices
on the `5 K4` side. It deletes no ideal-base edge. The added degrees are

```text
K3 side:  4,4,4
K4 side:  2^6,0^14.
```

Equivalently, the final defect-graph degree multiset is `6^3,5^6,3^14`.
The scope is the fixed-base extension slice, not every 23-vertex graph with
that degree multiset.

## Exact count and quotient

Every active `K4`-side vertex has two neighbors in the `K3`, so it is encoded
by the one `K3` vertex it misses. Each of the three missing colors occurs
twice. The labeled count is therefore

```text
C(20,6) * 6!/(2!^3) = 3,488,400.
```

The primary quotient enumerates all 3,488,400 labeled states. Inside each
`K4`, an `S4` orbit is completely determined by the three missing-color
counts. It sorts the five block-count triples for `S5` and minimizes over all
six missing-color permutations for `S3`. Equal canonical keys are therefore
equivalent under exactly

```text
S3 x (S4 wr S5), order 5,733,089,280.
```

This gives exactly 20 orbits. Their labeled-size histogram is:

```text
2160:1, 5760:1, 8640:3, 12960:1, 34560:1,
77760:1, 103680:3, 138240:1, 184320:1,
207360:3, 276480:1, 414720:1, 552960:1, 829440:1.
```

An independent quotient check enumerates 3,310 ordered block profiles,
whose multinomial multiplicities again sum to 3,488,400. Connected components
under adjacent block swaps and adjacent color swaps give the same 20
components and the same orbit-size multiset. Collapsing block order before
the color quotient independently gives 61 block multisets.

## Coverage relative to earlier exact searches

Exactly one orbit belongs to the original distinct-block-pair
three-`K2,2` connector family. Allowing repeated block pairs covers one more,
so the extended connector harness covers 2 of 20 and this audit adds 18
previously untested classes.

Exact maximum edge overlap with the published orbit gives the minimum
remove/add swap-radius histogram:

```text
radius 0: 1
radius 2: 3
radius 3: 3
radius 4: 8
radius 5: 4
radius 6: 1
```

Thus seven classes intersect the already-screened radius-at-most-three
family and thirteen do not. The union of connector-reuse and radius-at-most
three covers eight classes, leaving twelve outside both earlier families.
The independent Python replay recomputes the minimum radius from block
profiles and reproduces this histogram.

## Exact determinant screen

All 20 quotient representatives are positive definite. Four have determinant
strictly above the squared frontier. The largest is

```text
7,914,907,574,363,750,400,000,000,000,000,
```

and is not a square. There are exactly three square orbits:

```text
2,779,447,296,000,000^2   (one frontier tie)
2,739,929,088,000,000^2   (two distinct subfrontier orbits)
```

Every root is divisible by `2^22`, but no square is strictly above the
frontier. The strict factor-routing snapshot is therefore complete and empty;
the candidate was not changed.

The determinant engine uses four-prime centered CRT. The supported principal
minors are bounded by `727^(23/2) < 2^110`, while the CRT modulus exceeds
`2^123`, so reconstruction is exact and unique. A separate Python
arbitrary-precision Bareiss replay reconstructs all 20 Grams from their edge
sets, reproduces every determinant, checks every Sylvester minor, rebuilds the
20-class quotient, and matches the screen aggregates.

## Research-square gates

A separate research-only snapshot retains all three square orbits without
presenting the tie or subfrontier roots as strict improvements.

Exact Hasse--Minkowski screening finds no rational obstruction for any of the
three. The mod-3 sign-column-shell gate then:

- passes the known frontier orbit with shell size 1,382 and equal span ranks;
- exactly rejects the first subfrontier orbit with shell size 224 and ranks
  `62/63`;
- exactly rejects the second subfrontier orbit with shell size 396 and ranks
  `104/105`.

The retained shell report includes explicit mod-3 separating certificates.
Consequently neither subfrontier square can be the Gram matrix of a real
order-23 sign matrix.

## Reproduce

```sh
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wshadow -Werror \
  research/gram_published_degree_slice.cpp \
  -o build/research/gram_published_degree_slice

build/research/gram_published_degree_slice --self-test

build/research/gram_published_degree_slice \
  --output \
    runs/direct-search/gram-published-degree-slice-20260728/report.json \
  --route-snapshot \
    runs/direct-search/gram-published-degree-slice-20260728/route.json \
  --research-square-snapshot \
    runs/direct-search/gram-published-degree-slice-20260728/research-squares.json

python3 research/gram_published_degree_slice_replay.py \
  --report \
    runs/direct-search/gram-published-degree-slice-20260728/report.json \
  --output \
    runs/direct-search/gram-published-degree-slice-20260728/independent-bareiss-replay.json

python3 research/gram_hasse.py \
  --snapshot \
    runs/direct-search/gram-published-degree-slice-20260728/research-squares.json \
  --all-hits \
  --output \
    runs/direct-search/gram-published-degree-slice-20260728/research-squares-hasse.json

build/research/gram_shell_filter \
  --snapshot \
    runs/direct-search/gram-published-degree-slice-20260728/research-squares.json \
  --all-hits \
  --allow-subfrontier-research \
  --output \
    runs/direct-search/gram-published-degree-slice-20260728/research-squares-shell.json
```

The standard shell filter rejects roots at or below the frontier. The retained
research-square shell artifact was produced with its explicit opt-in
subfrontier-research mode; this does not alter strict promotion behavior.

## Frozen hashes

```text
enumerator source
48c1bc4b3fa680bb080be5af6b0606eb61176485e49009d1f1ebcd53ca48fdd5

independent replay source
d078751a772e2c144b6db6b4e448b3832adaef628d23ea152f04b2c3977937a2

optimized enumerator
b2789a39e30e9d4560aff13154af3380c32aa12459b6195d2bfc55beb96bc0de

primary report
a1a288d28c1431e861f2f2c4f12234f95c26f463c1be28c6df48927bce3490e7

strict route
17e10e2cc9c58cb080108661b4a37de56c93c15fac87127bc82c57a59e84fecd

research squares
608d3bf712fa78d68771f0bc22f532bb7ac91ec7fd6bdd5d8b3adc000cf65de4

independent Bareiss replay
6a354cbfe16641d818fbb0aede8b6bcdefba3e137bb7c16e776866281e7317e5

research Hasse report
421fe3afbc952c8caad6ee44484dec618ad30c48bf573965e4eb8929f434a4da

research shell report
772349ca4cadcc07cf423d7aa7a83dbc07aa0cc6a322668c0bea8f2d7356fcc7
```

The complete command and artifact manifest is
`runs/direct-search/gram-published-degree-slice-20260728/provenance.json`.

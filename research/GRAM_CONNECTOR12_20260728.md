# Exact three-connector Gram audit

**Run date:** 2026-07-28

This audit closes one structured part of the 12-edge gap between the ideal
order-23 Ehlich defect graph and the published frontier defect graph.

## Family

Start from

```text
K3 disjoint-union 5 K4
```

and choose three distinct unordered pairs of its six blocks. On each chosen
block pair, add one `K2,2`: choose two vertices in each endpoint block and add
all four cross edges. Connectors may share a block and may reuse vertices
inside that block. Distinct block pairs ensure that the result always has
exactly 12 distinct added edges.

There are exactly

```text
12,072,240
```

labeled configurations:

```text
[x^3] (1 + 18 x)^5 (1 + 36 x)^10
  = 12,072,240.
```

The five factors with weight 18 are the `K3`--`K4` block pairs
(`C(3,2) C(4,2) = 18`). The ten factors with weight 36 are the `K4`--`K4`
block pairs.

`gram_connector12.cpp` quotients this family exactly by the base graph
automorphism group

```text
S3 x (S4 wr S5), of order 5,733,089,280.
```

For each configuration it enumerates the six connector orders and eight
endpoint orientations. Active `K4` blocks receive first-occurrence labels,
and the ordered subset sequence at each block is minimized over all vertex
permutations of that block. The least resulting 42-bit word is a complete
canonical key: equality supplies explicit block, vertex, connector, and
endpoint bijections in both directions.

The 12,072,240 configurations reduce to exactly

```text
73 orbits.
```

The retained orbit-size histogram sums back to 12,072,240, and every orbit
size divides the full base automorphism group order.

## Published reproduction

The documented frontier construction uses three `K2,2` connectors from the
`K3` to three separate `K4` blocks, with the three different two-subsets of
the `K3`. Its canonical key is

```text
108190548787
```

and exact determinant is

```text
det(G) = 7,725,327,271,241,711,616,000,000,000,000
       = 2,779,447,296,000,000^2.
```

The published frontier factor independently passes `./arena verify`. Its
exact sign-column shell has size 1,382, contains every factor column, and has
no shell-span obstruction.

## Complete result

All 73 orbit representatives were evaluated with exact modular determinants.
Four have `det(G)` strictly above the squared frontier, and the maximum is

```text
7,862,828,089,422,643,200,000,000,000,000.
```

That maximum is not a square. Across the complete family:

```text
positive determinants             73
determinants above frontier^2       4
perfect-square determinants         1
frontier square ties                 1
perfect squares above frontier       0
qualified factor-routing survivors   0
```

Thus the published orbit is the only square Gram determinant in this
three-`K2,2`, distinct-block-pair family. No above-frontier candidate reaches
the Hasse, sign-shell, or factorization stages; their compatibility snapshots
are complete empty selections, not timeouts.

The determinant implementation uses four-prime CRT and symmetric
reconstruction. Every supported principal minor is bounded by
`727^(23/2) < 2^110` from its row norms, while the CRT modulus exceeds
`2^124`, making reconstruction exact and unique. An independent Python
Bareiss replay reproduced both the published square and the maximum
nonsquare determinant.

## Reproduce

```sh
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wshadow -Werror \
  research/gram_connector12.cpp -o /tmp/gram_connector12

/tmp/gram_connector12 \
  --output runs/direct-search/gram-connector12-complete-20260728.json \
  --route-snapshot \
    runs/direct-search/gram-connector12-complete-20260728-route.json

python3 research/gram_hasse.py \
  --snapshot \
    runs/direct-search/gram-connector12-complete-20260728-route.json \
  --all-hits \
  --output runs/direct-search/gram-connector12-complete-20260728-hasse.json

/tmp/gram_shell_filter \
  --snapshot \
    runs/direct-search/gram-connector12-complete-20260728-route.json \
  --all-hits \
  --output runs/direct-search/gram-connector12-complete-20260728-shell.json
```

The retained primary report SHA-256 is

```text
52ea3702c0f0ae8ce4bd715815d3d5d5be5980a61292179d4afc54084dd458de
```

Its complete run took 10.468 seconds on the audit machine.

# Exact targeted radius-4 frontier-basin closure

**Run date:** 2026-07-28

## Result

The frozen exact radius-3 catalog contains `1,958` `Aut(B0)` orbit
representatives, representing `148,736` labeled states, whose Gram
determinants are strictly above the squared frontier. For every such
representative, the screen enumerated all true outward moves:

```text
42 remaining base edges to delete
x 205 remaining base nonedges to add
= 8,610 radius-4 children per representative

1,958 x 8,610 = 16,858,380 parent-child transitions
```

These are transitions, not unique radius-4 states or orbits. The exact screen
found:

```text
positive-determinant transitions              16,858,380
determinants above the squared frontier          832,675
perfect-square transitions                           257
strict-above-frontier square transitions              29
distinct stored-coordinate states                       20
Aut(B0) route orbits                                      4
full graph/Gram permutation classes                       3
```

All four `Aut(B0)` route orbits were positive definite, had roots divisible by
`2^22`, and passed the exact rational Hasse screen. Their potential determinant
roots were:

```text
2,783,969,280,000,000   (two Aut(B0) orbits, one full class)
2,780,430,336,000,000
2,792,226,816,000,000
```

The exact sign-column shell/span filter rejected every class with a mod-3
certificate:

```text
route 0: shell 168, span 67, augmented rank 68
route 1: shell 168, span 67, augmented rank 68
route 2: shell   3, span  3, augmented rank  4
route 3: shell  10, span 10, augmented rank 11
```

Therefore none of these Grams can be the Gram matrix of a 23-by-23 sign
matrix. No matrix was constructed and the arena frontier remains
`2,779,447,296,000,000`.

The largest determinant encountered was the nonsquare positive-definite Gram
determinant

```text
8,109,553,344,734,822,400,000,000,000,000.
```

## Exact coverage boundary

Every radius-4 defect graph has 16 radius-3 parents: choose one of its four
deleted base edges and one of its four added base nonedges to undo.

Let a radius-4 state `C` have a radius-3 parent `P` with
`det(G(P)) > frontier^2`. The frozen radius-3 catalog contains one
representative `P0` of the `Aut(B0)` orbit of `P`. Some `g` in `Aut(B0)` maps
`P` to `P0`; the same `g` maps `C` to one of the `42 x 205` outward children
enumerated from `P0`. Thus this run covers exactly every radius-4
`Aut(B0)` orbit having at least one strict-above-threshold radius-3 parent.

It does **not** cover a radius-4 orbit whose entire set of 16 radius-3 parents
is at or below the threshold, and it is not a full radius-4 enumeration.

## Exact arithmetic and independent replay

The catalog is accepted only at frozen SHA-256

```text
b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed
```

The primary child screen applies each two-edge update to its radius-3 parent,
so the determinant lemma has endpoint dimension at most four. A separate full
replay applies every eight-edge radius-4 update directly to the reference Gram,
with endpoint dimension at most 16. Both use four proven 31-bit primes and
centered CRT. The row-norm bound

```text
727^(23/2) < 2^110
```

makes reconstruction unique under the greater-than-`2^123` CRT modulus.

The two implementations reproduced exactly:

- all transition statistics;
- every maximum and relevant key;
- all four determinant roots;
- byte-identical parent and route snapshots.

The screen also used full 23-by-23 elimination plus an independent fifth prime
on 129 deterministic transition samples, every relevant square, and the
maximum. GMP fraction-free Bareiss independently reproduced the maximum and
every relevant determinant. An independent Python Bareiss/Sylvester replay
confirmed all four stored route hits.

The four `Aut(B0)` route orbits collapse to three full defect-graph
isomorphism classes. Routes 0 and 1 are mapped by the vertex permutation that
swaps vertices 12 and 14 and fixes the other 21 vertices. Pairwise exhaustive
VF2 found the other two classes distinct. The certificate and pairwise matrix
are retained in `full-isomorphism-audit.json`.

## Reproduce

```sh
clang++ -std=c++20 -O3 -DNDEBUG \
  -Wall -Wextra -Wpedantic -Wshadow -Werror \
  $(pkg-config --cflags gmpxx) \
  research/gram_radius4_basin.cpp \
  -o build/research/gram_radius4_basin \
  $(pkg-config --libs gmpxx)

build/research/gram_radius4_basin --self-test

build/research/gram_radius4_basin \
  --screen \
  --catalog runs/direct-search/gram-radius3-orbits-20260728/catalog.tsv \
  --catalog-sha256 \
    b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed \
  --output \
    runs/direct-search/gram-radius4-basin-20260728/screen-parent-rank4.json \
  --parents-snapshot \
    runs/direct-search/gram-radius4-basin-20260728/parents.tsv \
  --route-snapshot \
    runs/direct-search/gram-radius4-basin-20260728/route-parent-rank4.json \
  --method parent-rank4 \
  --threads 8

build/research/gram_radius4_basin \
  --screen \
  --catalog runs/direct-search/gram-radius3-orbits-20260728/catalog.tsv \
  --catalog-sha256 \
    b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed \
  --output \
    runs/direct-search/gram-radius4-basin-20260728/screen-base-rank16.json \
  --parents-snapshot \
    runs/direct-search/gram-radius4-basin-20260728/parents-audit.tsv \
  --route-snapshot \
    runs/direct-search/gram-radius4-basin-20260728/route-base-rank16.json \
  --method base-rank16 \
  --threads 8

python3 research/gram_hasse.py \
  --snapshot \
    runs/direct-search/gram-radius4-basin-20260728/route-parent-rank4.json \
  --all-hits \
  --output runs/direct-search/gram-radius4-basin-20260728/hasse.json

build/research/gram_shell_filter \
  --snapshot \
    runs/direct-search/gram-radius4-basin-20260728/route-parent-rank4.json \
  --all-hits \
  --output \
    runs/direct-search/gram-radius4-basin-20260728/shell-parent-rank4-independent.json
```

## Frozen hashes

```text
source
a4a8d1829d88f33338d0e46fd71c0958997c9f56dc41cdcc6db1511dfaf0f351

included radius-3 trusted core
d4b0038d46c21b0a7af072a27cd3fdd5b3f29d6cc9ab7cae815d89f3ee60f512

optimized binary
87cb073989f8612f82ea424a5a5ea10d40994699cbd9eaf520598b7923cb331e

parent snapshots (both)
964f5c01cd3eb98dda5b90912afdff2116a5bb5a7300bab4c087839f68b6a9ae

primary screen
fe93070b3106a27259162257755410e7e25f5bf868b8e9543b7a5cb71664c64f

independent base-rank16 screen
333de2141eb872b5e3085731a288f192b9d31f6a429dc157f12ea9f10b5977e5

route snapshots (both)
1b2c4abda6d7f4ddcb14c4cdf233ae0d86835a69d4f76e89a37d2ff49220a3e6

Hasse report
75e2d203b33d1d444c8ecfbc6c430b2aa24b029672bfb124a07d218dc45bf377

shell report
882ff4eba9d44e33555c36833f2062ad8cd26594a29ea381c3fb00e9f2b60206

full-isomorphism audit
98a6e3401d2d613ec3c9b914b352af3ade34a92202a11a8bb17120761e18cfcc
```

# Structured order-23 searches: Hadamard descendants and order-22 borders

Date: 2026-07-28

This note audits the order-24 Hadamard coverage in the direct-search campaign
and records two finite, provenance-bound structured searches. All determinant
comparisons below are exact integers. No result in this note is an optimality
or inequivalence claim.

## What the original H24 generator covered

`research/hadamard24_seed.py` parses one of the 60 order-24 Hadamard
representatives, proves `HH^T = 24I`, and deletes one caller-selected row and
column. It emits one matrix per invocation. The main 60-class campaign invoked
it at deletion `(1,1)` for every class. Only three other pivots were retained
before this audit:

```text
class 14: row 7,  column 13
class 42: row 17, column 5
class 51: row 11, column 19
```

Thus, the original campaign covered all 60 Hadamard matrix equivalence-class
representatives but did not cover all designs obtained by dephasing those
matrices at different row/column pivots.

The local Spence catalog,
`runs/direct-search/reference-data/spence-hadamard.24`, has SHA-256
`cfdf2082c6a70d65a3e84b1167514cbe12a503a3e6d9cb0048b37f6c27fbd940`.
After dephasing each Mendeley representative at `(1,1)`, all 60 cores match
the 60 catalog cores exactly, including their ordering. The catalog's
per-class counts sum to 1,106 non-isomorphic Hadamard
`2-(23,11,5)` designs. Therefore the 60 default cores sampled at most
`60 / 1106` of this design family; the three one-off pivots added at most
three more design types.

The validation artifact is:

```text
runs/direct-search/h24-catalog-audit/summary.json
```

The classes with the most design descendants are useful diversity targets:

| H24 class | non-isomorphic designs |
| ---: | ---: |
| 35 | 144 |
| 19, 36, 38 | 72 each |
| 30, 39 | 48 each |
| 20, 21, 54 | 36 each |
| 8, 42 | 30 each |

For comparison, classes 14, 42, and 51 contribute 2, 30, and 2 designs.

Every direct order-23 minor of an order-24 Hadamard matrix has determinant
magnitude `24^11 = 1,521,681,143,169,024`. This follows from
`H^-1 = H^T / 24`: every cofactor has magnitude
`|det(H)| / 24 = 24^12 / 24`. Different pivots do not improve that starting
score, but they provide different factorizations and different trajectories
once local search leaves the Hadamard-minor plateau.

Explicit closed-quadruple switching at order 24 is not a larger seed universe:
a switched matrix is still in one of the 60 complete order-24 equivalence
classes. Enumerating every dephasing pivot of every representative covers the
associated 1,106 design descendants without needing to rediscover class
representatives by switching. Switching can still be useful as a cheaper way
to move between representations during a search.

## Multi-pivot Hadamard campaign

`research/h24_deletion_campaign.py` enumerates specified H24 classes and all
576 row/column pivots per class. It:

1. proves every source is Hadamard with the existing strict parser;
2. checks the dephased `(1,1)` core against the Spence catalog;
3. derives each 23-by-23 minor deterministically;
4. runs the exact-checkpointed reactive tabu engine with a derived seed;
5. independently recomputes each retained score with Python Bareiss
   elimination;
6. atomically checkpoints a bounded elite set and a SHA-256-bound summary.

Representative command:

```sh
python3 research/h24_deletion_campaign.py \
  --classes 14,42,51 \
  --seconds-per-minor 0.25 \
  --seed-base 32000 \
  --top-count 30 \
  --output-dir runs/direct-search/h24-deletion-14-42-51
```

The first campaign completed:

```text
classes:       14, 42, 51
pivots:        1,728 / 1,728
search time:   432.0 aggregate seconds
wall time:     676.308 seconds
iterations:    130,881,009
exact checks:  2,355,467
best class:    42
best pivot:    row 13, column 12
best |det|:    2,694,905,856,000,000
frontier gap:  84,541,440,000,000
```

Best arena receipt:

```text
raw SHA-256:
a1bf5bea2e876af60b6fab424f5f03688d47bdbcd99046eba25b0186397d5d3d

normalized SHA-256:
4434926d8432081d411c67fd7680822a8fc932bcf9000fef75bf158535c27cc9

receipt SHA-256:
b186cee2e94df3c6b0eed1238b9d3808425d2224e0bfae3517fb499be549b5d7
```

The matrix is:

```text
runs/direct-search/h24-deletion-14-42-51/elites/class-42-r13-c12.matrix.txt
```

Its score ties a previously retained class-51-derived checkpoint but its raw
and sign-normalized hashes differ. Those hash differences do not establish
Hadamard inequivalence, so no novelty claim is made.

Two selected continuations added 26,167,245 moves and 408,888 exact checks
without improving their starts. The stronger continuation ran 144.314
seconds from the campaign leader and retained the same score.

The high-design-diversity campaign then covered all pivots of classes
19, 30, 35, 36, 38, and 39. Those six classes contain 456 of the catalog's
1,106 non-isomorphic designs:

```text
pivots:        3,456 / 3,456
search time:   345.6 aggregate seconds
wall time:     392.502 seconds
iterations:    98,988,517
exact checks:  2,119,448
best class:    36
best pivot:    row 8, column 3
best |det|:    2,580,524,128,272,384
```

A 180.009-second continuation from that shallow leader reached
`2,694,905,856,000,000` after 13,920,016 moves and finished after
36,415,919 moves and 569,014 exact checks without further improvement.
It independently reproduces the strongest score of the first campaign from a
different H24 class, but remains below the comparison frontier.

```text
runs/direct-search/h24-deletion-high-diversity/
  followup-class36-r8-c3.matrix.txt

raw SHA-256:
63038a1bc0a25acb4f16136b3973e9b059314d84079d1bd007c5b5497630ca31

normalized SHA-256:
715391f6b89de3eb89de8204519064fe04dc35c5a0cb60540f58f70fdee51276

receipt SHA-256:
001eb313f7adbb3088386a495489e04aebb9c65c2981267fabe53a32ceb86f28
```

Together, the two multi-pivot screens covered 5,184 pivots from H24 classes
containing 490 catalog design types. They performed 229,869,526 reactive
iterations and 4,474,915 independent exact checks. Their three continuations
added 62,583,164 moves and 977,902 exact checks. No strict frontier
improvement was found.

A final background campaign completed every pivot in the other 51 classes:

```text
runs/direct-search/h24-deletion-remaining-classes/summary.json
```

It covered all 29,376 selected pivots whose classes contain the remaining 616
catalog design types, using 330,843,073 reactive iterations and 10,178,814
exact checks at 0.05 seconds per pivot. It found an exact frontier tie from
H24 class 9 with row 10 and column 24 deleted, and no strict improvement.
Together the three campaigns cover all 34,560 row/column pivots and all 1,106
catalog design descendants represented by the pinned 60-class catalog.

## Exhaustive bordering of the order-22 corpus maximum

This is structurally independent of the H24-minor search. The extractor:

```sh
python3 research/corpus_matrix.py \
  --order 22 \
  --output runs/direct-search/order22-border/mendeley-order22.matrix.txt
```

reads the sole order-22 matrix from the local Mendeley maximal-determinant
corpus. Provenance:

```text
corpus SHA-256:
4e963528533cd0a839bb9073729095e4274c2747416a205a9b83935ce9dd3c8a

extracted core SHA-256:
c7684c22aa4dc37fa18617e186c292268d2cb5f70b2a04268dd1c4b008c04f93

|det(core)|:
409,600,000,000,000
```

For a fixed core `B`, write a border as

```text
A = [ B  x ]
    [ y' c ]
```

Then, exactly,

```text
det(A) = det(B)c - y' adj(B)x.
```

Fixing the first sign of `x` removes the global sign redundancy. For each of
the resulting `2^21` columns, the optimal row and corner are immediate:
choose `c = sign(det(B))` and choose each `y_i` opposite in sign to
`(adj(B)x)_i`. The exact objective is therefore:

```text
|det(B)| + sum_i |(adj(B)x)_i|.
```

`research/order22_border.cpp` computes the adjugate with integer cofactors,
checks `B adj(B) = det(B) I`, traverses all `2^21` columns in Gray-code
order, and Bareiss-checks every promoted matrix.

```sh
c++ -std=c++20 -O3 -march=native \
  -Wall -Wextra -Werror -pedantic \
  research/order22_border.cpp \
  -o build/research/order22_border

build/research/order22_border \
  --start runs/direct-search/order22-border/mendeley-order22.matrix.txt \
  --output runs/direct-search/order22-border/best.matrix.txt \
  --log runs/direct-search/order22-border/screen.jsonl
```

The exhaustive result was:

```text
assignments:  2,097,152 / 2,097,152
best |det|:   2,465,792,000,000,000
best columns: 10 up to global sign
elapsed:      0.051 seconds
```

Arena verification:

```text
raw SHA-256:
5a70b81299eacc74176f2b7ba9e0eda0c1f2bd83d63850b1057419032114d8c2

normalized SHA-256:
7366f2997cb82af6d79459eea2e61c7cd63c4e005fa7fcceb65ad8d145698cfe

receipt SHA-256:
3fc9b4f99168a9553ca760b23b8ab556f8a11215d74feca7206ae0e58b38f745
```

An ASan/UBSan build reproduced the same matrix. This closes the complete
one-row/one-column border family around this specific order-22 core; it does
not cover the second known order-22 information-matrix class or perturbed
order-22 cores.

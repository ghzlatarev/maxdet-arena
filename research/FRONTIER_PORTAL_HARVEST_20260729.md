# Frontier portal harvest and exact shell closure (2026-07-29)

## Outcome

The bounded fixed-Gram harvester pilot produced no ninth local H/HT class and
no new Gram class. It retained two exact frontier ties:

```text
|det| = 2,779,447,296,000,000
```

Both ties passed exact Gram reconstruction, integer Bareiss, `./arena verify`,
and the pinned `pynauty==2.8.8.1` H/HT/Gram classifier. They reduce to already
known local H/HT classes `b64c3309...` and `eb138a06...`.

A separate polar/Douglas–Rachford campaign completed later and did find two
additional local H/HT classes. Exact alignment places both in this same
published Gram and the surviving `(0,2,4)` triple orbit. That later result is
consistent with the theorem below: the theorem eliminates the other three
triple orbits, but does not classify or exhaust factors in `(0,2,4)`.

The more important result is exact: the three previously unresolved
size-six-shell triple orbits are impossible for this fixed Gram. Two sparse
integer identities exclude their representatives, and in fact exclude 14 of
the 20 size-three subsets without invoking an LP or SAT solver.

## Claim boundary

The fixed-Gram pilot's frozen input baseline was the explicit local corpus and
pinned classifier state:

```text
11 H classes / 8 H/HT classes / 1 normalized row-Gram class
```

The later polar/Douglas–Rachford campaign expanded that local corpus to
`14 H / 10 H/HT / 1 normalized row-Gram class`. Neither count is a literature
classification or priority claim. The exact result here closes three
column-factor slices of one fixed row Gram. It does not prove that the ten
currently local H/HT classes exhaust that Gram's remaining factor slice, that
other frontier Grams do not exist, or that the current determinant is globally
optimal.

The sparse identities were discovered locally from a numerical separation
signal and then independently reduced to and checked as integer equalities.
No claim is made that the identities are new to the literature.

## Fixed-Gram portal harvester

`frontier_portal_harvest.py` searches the complete 1,382-column normalized
sign shell of the published frontier Gram. Exact orbit equations force factor
incidences

```text
3, 6, 6, 8
```

in shell orbits of sizes `6, 432, 432, 512`. A nonsingular normalized sign
factor cannot repeat a column: repeated normalized columns would be equal,
making the factor singular. The size-six incidence is therefore a three-mask
subset.

All eight then-known aligned local H/HT representatives select the same
canonical small-orbit triple `(0,2,4)`. The two classes found later by the
polar/Douglas–Rachford campaign align to that triple as well. The harvester:

- assigns deterministic global trial indices under `--shard-index` and
  `--shard-count`;
- prioritizes the other three canonical triple types before known-type
  controls;
- uses seeded positive HiGHS objectives that penalize frequently used masks or
  overlap with a rotating known support;
- deliberately adds no exact raw-support no-goods;
- retains only exact Gram/Bareiss frontier factors and gives every retained
  factor an arena receipt and pinned H/HT/Gram classification.

The final provenance-bound pilot completed eight trials. Each missing triple
type `(0,1,2)`, `(0,1,3)`, and `(0,2,5)` received a rarity and an anchor
objective. Those six models returned HiGHS infeasible, while both `(0,2,4)`
controls returned exact frontier factors. The solver statuses alone were
treated only as bounded floating-MILP evidence; the next section supplies the
independent proof.

The two control supports had minimum symmetric difference `36` from the eight
aligned representatives but still collapsed to known H/HT classes. This is
direct evidence that raw aligned-support distance is a weak novelty proxy and
explains why raw-subset no-goods previously reduced throughput without
reliably excluding equivalence classes.

Authoritative pilot manifest:

```text
runs/direct-search/frontier-portal-harvest-20260729/pilot-wave2/manifest.json
campaign receipt:
676451d6c691583a926bcc39b30ccb557875b943f2d145aae23a7d68b4277ffb
```

The manifest records the actual subprocess environment:

```text
Python 3.9.6 / NumPy 1.22.4 / SciPy 1.13.1 / pynauty 2.8.8.1
```

## Exact integer certificate

Let `s` be a normalized shell sign column, with zero-based coordinates. Define

```text
L1(s) = -(s[9] + s[10]) (s[15] + s[16] + s[17] + s[18])
L2(s) =  (s[3] + s[4])  (s[15] + s[16] + s[17] + s[18]).
```

Each expression is an eight-term integer linear combination of the
off-diagonal outer-product coordinates `s[r]s[c]`. Direct enumeration of the
hashed complete shell gives

```text
L1(s) = L2(s) = 0
```

for every one of the 1,376 columns outside the size-six orbit. Summing the
same eight Gram entries in the target gives the forced factor totals

```text
sum L1 =  8
sum L2 = -8.
```

Only the three selected size-six columns can contribute. The missing
representatives contradict these totals:

| fixed triple | functional | fixed contribution | target | residual required from the 1,376-column kernel |
| --- | --- | ---: | ---: | ---: |
| `(0,1,2)` | `L1` | `24` | `8` | `-16` |
| `(0,1,3)` | `L2` | `8` | `-8` | `-16` |
| `(0,2,5)` | `L1` | `24` | `8` | `-16` |

The remaining columns contribute exactly zero under the relevant functional,
regardless of their coefficients. They therefore cannot supply `-16`. This is
an integer linear inconsistency, stronger than an integrality or nonnegativity
obstruction.

Checking all `C(6,3)=20` small-orbit subsets leaves exactly

```text
(0,1,4), (0,2,4), (0,4,5), (1,2,3), (1,3,5), (2,3,5).
```

These are the six orientations in the observed canonical orbit; the other 14
subsets are rejected by `L1`, `L2`, or both.

The standalone checker uses Python arbitrary-precision integers only:

```sh
python3 research/frontier_portal_farkas_check.py --stdout
```

Pinned artifacts:

```text
checker:
research/frontier_portal_farkas_check.py
SHA-256 5bbb7efc1f95083449cea7995f81dec154914596cbbd9f15d7d1d8f991346c3d

report:
runs/direct-search/frontier-portal-harvest-20260729/exact-farkas-report.json
SHA-256 4cc981632f575974e59475301fbe218df5165ab1f67529ffbb52ae23c38855e2
```

The report binds the complete shell and the independent exact `3,6,6,8`
orbit-count certificate by SHA-256. Regenerating it with `--stdout` is
byte-identical to the retained report.

## Reproduction and next use

Run a fresh deterministic shard with a distinct output directory:

```sh
python3 research/frontier_portal_harvest.py \
  --output-dir runs/direct-search/frontier-portal-harvest-NEW/shard-0 \
  --shard-count 8 --shard-index 0 \
  --trials-per-shard 8 --time-limit 90 --jobs 2 \
  --classifier-python /path/to/pynauty-2.8.8.1/bin/python
```

The exact closure means further fixed-Gram factor harvesting should search only
the surviving triple orbit. The later polar/Douglas–Rachford campaign already
revealed the ninth and tenth local H/HT classes there, and further classes may
remain; support distance alone should not be used as their novelty certificate.
To find a new Gram portal, the search must leave this fixed-Gram model; every
factor constructed here has the published row Gram by definition.

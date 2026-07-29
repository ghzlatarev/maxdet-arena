# Frontier factor-class expansion (2026-07-28)

## Outcome

The exact frontier did not move:

```text
|det| = 2,779,447,296,000,000
```

The useful result is factor diversity. Against the repo-wide augmented local
baseline produced before this campaign (`9` H-classes, `7` H/HT-classes), the
campaign retained:

1. one missing transpose partner, raising the corpus to `10 H / 7 H/HT`;
2. one genuinely new local H/HT class from shell MILP seed `410503`, raising
   the final corpus to **`11 H / 8 H/HT`**.

“New” here means absent from the explicit local corpus under the pinned
classifier. It is not a literature-novelty claim. All `27` matrices in the
final audit remain in one normalized row-Gram class.

## Retained direct-search seeds

| class | matrix | internal arena receipt SHA-256 |
| --- | --- | --- |
| H-only class `df0b9405...`; transpose H-partner `5f3d7a03...` | `runs/direct-search/frontier-factor-class-expansion-20260728/seeds/h-df0b940533f84c9d61ec8df73000b4fe3a646f073a6fd32f3826d72331ebedc0.matrix.txt` | `4a477cfa73413596460c27c575e02a6b5947682e8b7efd40cb09b571825caa9a` |
| new self-dual H and H/HT class `eb138a06...` | `runs/direct-search/frontier-factor-class-expansion-20260728/seeds/h-eb138a06ec638735c34bdacf77bd1cdd869c5d2fbc3450be25d63af0cde1a134.matrix.txt` | `f7a77f60352d225aefa17e79651d548e18f34c86deab0513cab7d4bd2904f6d7` |

Both independently pass `./arena verify`. Duplicate generated factors were
classified and discarded; their exact solver metadata remains in the campaign
directory.

## Exact shell model and symmetry

The published frontier Gram has an exact normalized sign-column shell of
size `1,382`. Factor multiplicities are Boolean: two equal normalized columns
would make a 23-by-23 factor singular, contradicting the positive Gram
determinant.

The exact Gram automorphism group has order `442,368`, with `15` generators.
Its shell action has orbit sizes:

```text
6, 432, 432, 512
```

The action on unordered row pairs has `16` orbits. Summing the `253` Gram
equations over those pair orbits, and adding the 23-column equation, gives a
17-by-4 integer system of rank four. Exact elimination uniquely forces factor
incidences:

```text
3, 6, 6, 8
```

The twenty triples in the size-six orbit split into four automorphism orbits,
represented by indices:

```text
(0,1,2), (0,1,3), (0,2,4), (0,2,5)
```

The current exact derivation, dependency versions, source SHA, retained
calibration factor, and arena receipt are bound by:

```text
runs/direct-search/frontier-factor-class-expansion-20260728/calibration/orbit-proof-v4.json
```

## Bounded searches

### Orbit-fixed CP-SAT

Twelve substantive trials used exact Boolean constraints, forced orbit counts,
fixed small-orbit triples, and then class-hinted overlap diversity:

```text
verified factors     2
new H / H/HT         0 / 0
UNKNOWN timeouts    10
INFEASIBLE proofs    0
summed solver time   1,226.275324 s
```

The known `(0,2,4)` type calibrated in `1.520 s`. The other three triple types
were each `UNKNOWN` after 120 seconds. An overlap-at-most-16 run found a known
`4072ead4...` factor in `24.664 s`; the remaining overlap and aligned-hint
trials timed out. None of those timeouts proves an empty slice.

### Feasible-first HiGHS MILP

The successful formulation used:

```text
1,382 Boolean variables
strictly positive seeded objective
mip_rel_gap = 1.0
exact outer-product and Bareiss checks after a feasible solve
```

The retained claim is feasibility plus exact verification. HiGHS status 0 and
a nonzero reported MIP gap are not interpreted as objective optimality.

Across twelve 45-second seeds in three parallel waves:

```text
verified factors      5
new H / H/HT          1 / 1
UNKNOWN timeouts      7
INFEASIBLE proofs     0
process-real total  414.06 s
parallel wall       137.37 s
```

The five H-certificates were `b584c923...`, `4072ead4...`, `eb138a06...`,
`ff1b5d37...`, and `9035bdf2...`. The third is the retained new class.
A four-seed exact raw-subset no-good wave returned no factor; removing those
cuts restored feasibility throughput. The last four-seed wave yielded two
known classes and no new class, so the stream was stopped rather than
accumulating raw duplicates.

Seed `410503` was replayed byte-for-byte with the final pinned tool source:

```text
runs/direct-search/frontier-factor-class-expansion-20260728/calibration/milp-seed410503-replay-v2.json
```

## Verification and handoff

The authoritative compact report is:

```text
runs/direct-search/frontier-factor-class-expansion-20260728/campaign-report.json
```

The final pinned H/HT audit reports:

```text
input files       27
unique matrices   27
H-classes         11
H/HT-classes       8
Gram classes       1
```

at:

```text
runs/direct-search/frontier-factor-class-expansion-20260728/final-h-equivalence-audit.json
```

The `eb138a06...` factor is the higher-value new search seed because it adds a
new class even when transpose is allowed. It subsequently received two
completed core-adjugate arms and an exhaustive `2^32` connector search; see
`EB138A_SEARCH_20260728.md`. The `df0b9405...` transpose partner remains
available as an orientation-control seed.

The candidate was not changed because no strict score improvement appeared.

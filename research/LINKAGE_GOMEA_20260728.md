# Exact linkage-GOMEA pilot (2026-07-28)

## Result

The corrected five-minute pilot did **not** beat the configured order-23
frontier:

```
2779447296000000
```

The retained matrix is an already-known normalized frontier representative,
not a novelty claim. `./arena verify` reports:

- exact determinant: `-2779447296000000`
- raw/normalized SHA-256:
  `efea293d313742c5e0b72167d59bccb5f25384fd51c87dd65e2634e95a6b40f7`
- receipt SHA-256:
  `5ff4b239e7c9a1c1704eec5b21233d8b8d38a481bec3dcd669046aab64dff26a`

All 99 retained archive matrices also passed `./arena verify`.

## Method

`research/core_linkage_gomea.cpp` is a standalone exact
gene-pool-optimal-mixing pilot over the dephased 22-by-22 binary core.
It:

- seeds the 24 aligned neutral-network frontier matrices and one separately
  aligned H2/QUBO frontier matrix;
- learns a normalized-mutual-information linkage tree every four
  generations;
- includes core rows, core columns, all known neutral generator supports,
  and the transported reference-to-H2 support in the linkage family;
- copies donor subsets and accepts only exact nondecreasing Bareiss
  determinants;
- applies exact adjugate-guided one-bit polishing;
- rejects offspring or polished results that duplicate another live lineage;
- uses exact nondecreasing archive resets to preserve live diversity; and
- atomically checkpoints the best matrix, JSONL log, final JSON summary, and
  exact archive.

Objective checks are differential-tested against an independent determinant
modulo `1000000007` and the exact 23-by-23 sign/core determinant quotient.

## Coordinate alignment and bridge

Raw dephasing does not align row/column labels across H-classes. The QUBO
tie was therefore aligned to
`tie-000-0000-3fa9fd308dc5.matrix.txt` with signed row/column permutations:

- raw Hamming distance before alignment: 263
- aligned Hamming distance: 75
- exact determinant preserved: `2779447296000000`
- aligned matrix SHA-256:
  `854e37a0f682f3f950259b20120e3bcee315fc88a88178233d34d62b5998d2e8`

The same monomial transform was applied to the Orrick reference. The
transported reference and aligned H2 seed differ in exactly 12 sign entries,
and—importantly—in exactly 12 dephased core variables too. Both endpoints
pass `./arena verify`.

Artifacts:

- `runs/linkage-gomea-alignment-20260728/alignment.json`
- `runs/linkage-gomea-alignment-20260728/bridge.json`
- `runs/linkage-gomea-alignment-20260728/seventh-aligned.matrix.txt`
- `runs/linkage-gomea-alignment-20260728/reference-transported.matrix.txt`

## Harness quirks found and fixed

1. The twelve known neutral sign-matrix generators do not all remain
   12-variable moves after dephasing. Six induce 12-bit core supports; the
   six touching the gauge border induce 31-bit core supports. The solver now
   derives and strictly validates the `6 + 6` split instead of dropping half
   the known structure.

2. Merely validating the transported reference did not put the opposite
   bridge endpoint in the donor pool. The corrected solver adds it as an
   archive-only frontier donor, preserving the requested 25 population seeds,
   and chooses the exact opposite endpoint for the bridge subset.

3. Unrestricted equal-score drift collapsed the first pilot from 32 live
   cores to 2. The corrected solver rejects equal, improving, and polished
   duplicate outcomes and performs nondecreasing diversity resets from the
   archive. The second pilot stayed at 32/32 unique live cores.

## Baseline negative calibration

`runs/linkage-gomea-pilot-20260728-seed33007/`

- elapsed: 300.034735 s
- exact evaluations: 5,799,978
- generations: 576
- accepted improvements / ties: 58 / 5,640
- bridge changed attempts / accepts: 29,899 / 0
- final population uniqueness: 2 / 32
- archive: 85
- promotions: 0

This run is retained as the negative control demonstrating collapse and the
missing bridge endpoint.

## Corrected five-minute pilot

`runs/linkage-gomea-diverse-pilot-20260728-seed33013/`

- elapsed: 300.030419 s
- exact evaluations: 4,207,670
- generations: 604
- accepted improvements / ties: 63 / 66
- transported bridge accepts: 7
- rejected duplicate equal / improving / polished outcomes:
  2,341 / 3,620 / 891
- exact nondecreasing diversity resets: 2,398
- final population uniqueness: 32 / 32
- archive: 99
- differential checks: 92
- promotions: 0

Primary artifact hashes:

- best matrix:
  `efea293d313742c5e0b72167d59bccb5f25384fd51c87dd65e2634e95a6b40f7`
- summary:
  `9f121b2013a9d59b277b13570025291d003e5c67fb2d10db0ed5f587114a383f`
- JSONL log:
  `590f6c0824c280f28db684ff984904b7131cf3943d28796cbd022fc6b8daf80d`
- archive manifest:
  `5d00819d9b37f6f587b9f4767aa3ad3de3c69f1ec105779491ae13cb03bf6c5c`

## Validation

The source compiled with:

```
clang++ -std=c++20 -O3 -march=native -Wall -Wextra -Wshadow -Werror -pedantic
```

It also passed an AddressSanitizer + UndefinedBehaviorSanitizer run, strict
NaN/input rejection, rejection of an unaligned bridge whose 12 raw flips
induce 88 core flips, exact random donor-subset differential checks, JSON
parsing, fresh-output enforcement, and arena verification.

No H/HT novelty audit was triggered: there was no score promotion and this
report makes no archive-equivalence or novelty claim.

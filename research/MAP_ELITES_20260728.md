# Exact quality-diversity pilot (2026-07-28)

## Outcome

The corrected five-minute pilot did not beat the order-23 frontier:

```text
2779447296000000
```

The retained best independently passed `./arena verify` with normalized
SHA-256
`efea293d313742c5e0b72167d59bccb5f25384fd51c87dd65e2634e95a6b40f7`
and receipt
`5ff4b239e7c9a1c1704eec5b21233d8b8d38a481bec3dcd669046aab64dff26a`.

## Method

`research/core_map_elites.cpp` combines exact extremal-optimization repair
bursts with a MAP-Elites-style archive. A niche is keyed by the nearest of
three supplied frontier seeds, binned Hamming distance, the number of nearly
balanced rows and columns, and binned off-diagonal row-Gram energy.

Parent tournaments alternate between underused niches and high-determinant
niches. Three quarters of parents receive an exact 4-to-24-bit destroy step;
the best of twelve nonsingular destroy candidates is observed and archived
before a 512-move EO repair burst. Every accepted state update uses the
differential-tested exact adjugate kernel, and archive quality is the exact
integer determinant.

The motivation is the quality-diversity idea that structurally different
high-quality solutions can serve as stepping stones rather than allowing one
lineage to dominate. This is supported in other combinatorial settings by
[Qian, Xue, and Wang](https://arxiv.org/abs/2401.10539). The repair dynamics
adapt the rank-based downhill exploration of
[Boettcher and Percus](https://arxiv.org/abs/cond-mat/0010337).

## Harness bug found before the final run

The first pilot, seed `35007`, was stopped at 59.9 seconds after an independent
audit found that the selected exact destroy endpoint was not checked for a
frontier tie or strict promotion before the first EO repair move. A winning
kick could therefore have been moved away and lost.

The invalidated artifact is explicitly marked:

```text
runs/map-elites-pilot-20260728-seed35007/INVALIDATED.md
```

The corrected engine observes, checkpoints if better, and archives the exact
destroy endpoint before any repair move. Seed `35007` is excluded from every
campaign total below.

## Corrected pilot

Artifact:
`runs/map-elites-pilot-20260728-seed35011/`.

- elapsed after startup/self-test, including archive export: `300.011643 s`;
- EO moves: `7,590,400`;
- exact repeated direction evaluations: `3,905,433,389`;
- exact destroy candidates: `133,068`;
- occupied structural niches: `1,933`;
- archive replacements: `3,161`;
- exact state rebuilds: `14,826`;
- exact determinant / adjugate-identity / hash checks:
  `926 / 26,841 / 926`;
- raw dephased frontier cores seen on the trajectory: `37`;
- retained frontier elites among the top 32 outputs: `18`; and
- strict promotions: `0`.

All 32 exported archive matrices independently passed `./arena verify`.
Direction counts repeat neighborhoods and are not unique-matrix counts.
The 37 trajectory cores are exact binary-core identities, but only the 18
frontier elites selected by the one-per-niche archive were exported.

## Equivalence audit

The pinned `pynauty==2.8.8.1` classifier compared the 18 retained frontier
elites with the frozen local corpus and the previously found QUBO class.
Relative to that seven-class augmented corpus:

- H-equivalence classes increased from 7 to 9;
- H/HT classes, allowing transpose, stayed at 7; and
- normalized row-Gram classes stayed at 1.

Thus the two additional H-only classes are transpose partners of already
known local classes. They are useful evidence that the diversity archive
crossed orientation basins, but they are not new MaxDet equivalence classes.
The compact audit is
`runs/map-elites-pilot-20260728-seed35011/h-equivalence-summary.json`.

## Validation and provenance

The engine passed a strict warning build, exact mixed-move differential tests,
an aggressive all-kick/archive-every-move smoke test, and an
AddressSanitizer/UndefinedBehaviorSanitizer run. The retained best and all 32
archive exports passed the trusted arena verifier.

Hashes of the exact executed artifacts:

```text
78c4d4ff82cac2e6c0d782ca87408cdebab5c48d06512bdcd22856c847c20afb  research/core_map_elites.cpp
0bbfc443d1f971a020481b36df00e0a7a2daaeddf51e25b969d3c5d87b52e3e0  research/core_extremal_optimization.cpp
182d19b0f76d3c1f15c816a0745c2684196807b10afc6b14041f39940353ee0e  build/research/core_map_elites
38b989ffb81f2972488d4d375dabce5faa832e1a268828936450093250d6bd2c  summary.json
d1f442a5136499677999087e3ec1f87b98add47a6862a46c6f28023085660f35  events.jsonl
```

Representative command:

```sh
build/research/core_map_elites \
  --seed-matrix runs/direct-search/neutral-cycle/two-cycle-union-29952/tie-000-0000-3fa9fd308dc5.matrix.txt \
  --seed-matrix runs/direct-search/neutral-cycle/two-cycle-union-29952/tie-001-0001-edbf3e96a40f.matrix.txt \
  --seed-matrix runs/qubo-trust-pilot-20260728-seed31003/best-proposal.matrix.txt \
  --output runs/map-elites-pilot-20260728-seed35011/best.matrix.txt \
  --archive-dir runs/map-elites-pilot-20260728-seed35011/archive \
  --log runs/map-elites-pilot-20260728-seed35011/events.jsonl \
  --summary runs/map-elites-pilot-20260728-seed35011/summary.json \
  --seconds 300 --heartbeat 15 --burst-moves 512 \
  --archive-stride 64 --differential-rounds 24 --seed 35011
```

This is a stochastic negative search result. Archive niches are not
equivalence classes, the local classifier is not a literature survey, and the
pilot makes no optimality or literature-novelty claim.

# QD-selected exact cubes and disjoint connector

Date: 2026-07-28

## Outcome

No strict improvement was found. The verified order-23 frontier remains:

```text
2,779,447,296,000,000
```

This tranche added `12,884,901,888` exact assignment-visits:

| Campaign | Logical cubes | Evaluator leaves | Assignment-visits | Outcome |
| --- | ---: | ---: | ---: | --- |
| QD-selected core supports | 48 × 27-bit + 1 × 31-bit | 64 | 8,589,934,592 | no gain |
| Elite 004↔aligned-021 connector | 1 × 32-bit | 32 | 4,294,967,296 | no gain |
| **New tranche** | **50** | **96** | **12,884,901,888** | **no gain** |
| **All FPM work to date** | **362** | **470** | **63,082,332,160** | **no gain** |

The total is an assignment-visit count. Different affine cubes can overlap.
The complete 31- and 32-bit cubes are internally disjoint partitions.

## QD-selected campaign

The hardened planner selected 12 distance/orientation-diverse exact-frontier
outputs from the corrected QD archive. The first five were independently
audited representatives of the five observed H/transpose classes: archive
ranks 003, 006, 011, 016, and 017.

For every center it built four distinct core-only 27-entry supports:

- exact one-flip gain;
- exact pair score;
- exact pair synergy; and
- a balanced single/synergy hybrid.

All support coordinates were restricted to the dephased 22×22 core, mapped to
one-based full-matrix rows and columns 2–23. Every support had 27 unique
coordinates and full GF(2) direction rank.

One of the two raw 31-entry frontier bridges was retained as a control.
The two candidates shared the same support, so only the elite 009↔010 bridge
was run. Its 16 fixed-prefix leaves exactly partitioned all `2^31`
assignments.

Final audit:

```text
plan SHA-256:                         7f888a20bbd9af254d40a75e702f8a3f55099a6ac6441b686368a87ce3ee7774
scheduled/evaluated tasks:            64 / 64
unrun tasks:                           0
assignment-visits:                     8,589,934,592
independent best replays:              64
independent top-K replays:             2,048
independent returned-tie replays:      1
truncated tie reports:                 0
strict promotions:                     0
```

The only logical nonzero frontier occurrence was the expected opposite
endpoint of the 31-bit bridge. It was not a new frontier matrix.

Evidence:

```text
runs/qd-selected-cubes-20260728-seed36002/plan.json
runs/qd-selected-cubes-20260728-seed36002/aggregate-report.json
research/qd_cube_campaign.py
```

The earlier `seed36001` directory is explicitly marked superseded. It was a
planning-only artifact; no evaluator task was run from it.

## Exact 32-bit connector

The strongest pair-selection result joined:

- frontier elite 004, score `2,779,447,296,000,000`; and
- subfrontier elite 021, score `2,726,756,352,000,000`.

A deterministic full 23×23 monomial alignment reduced their raw Hamming
distance from 289 to exactly 32. This is a full-matrix aligned connector, not
a 32-step transition observed by QD and not a dephased-core distance.

The 32 differing entries reconstruct the aligned endpoint exactly. Their
affine fingerprint is:

```text
91e8aa0c4f1ab6ddbbc0ea509620c7a14fb264299cc459631f8a87ccc1bce8d1
```

An exact dephased GF(2) audit gave rank 32/32 and an empty intersection with
each of the 312 earlier top-level affine cubes (200 initial, 110 LNPS, and
two deep-32 cubes). Therefore all `2^32 = 4,294,967,296` dephased states in
this connector are new relative to that pinned prior-cube corpus. This does
not claim disjointness from every order-23 matrix ever searched by other
methods.

The connector was split into 32 disjoint 27-bit leaves:

```text
completed leaves:                      32 / 32
assignment-visits:                     4,294,967,296
unique leaf fingerprints:              32
rerooted leaves:                        0
truncated tie reports:                  0
best score:                             2,779,447,296,000,000
strict promotions:                     0
```

The sole frontier assignment was global mask zero, the starting elite. The
aligned subfrontier endpoint was global mask `2^32-1` and ranked first among
the strictly subfrontier matrices.

Evidence:

```text
runs/direct-search/fast-principal-cube/
  qd-aligned-connector32-elite004-elite021-20260728-inputs/alignment.json
  qd-aligned-connector32-elite004-elite021-20260728/manifest.json
  qd-aligned-connector32-elite004-elite021-20260728/aggregate-report.json
  qd-aligned-connector32-elite004-elite021-20260728/affine-gf2-audit.json
  qd-aligned-connector32-elite004-elite021-20260728/provenance.json
  qd-aligned-connector32-elite004-elite021-20260728/provenance.sha256
research/affine_gf2_audit.py
research/fast_cube_connector_provenance.py
```

## Harness issues found and fixed

1. QD coordinates address the dephased core. Direct support generation now
   rejects first-row or first-column entries.
2. Equal affine fingerprints prove cube equality, not disjointness. The new
   GF(2) auditor detects exact intersections and containment after dephasing.
3. Every plan pins the driver, evaluator source/binary, imported helpers,
   arena verifier, QD summary, and prior-fingerprint evidence.
4. Shards use unique IDs and fresh attempt directories. Resume validates the
   immutable plan and every completed artifact before skipping work.
5. Every retained best, top-K record, and returned tie is reconstructed from
   its mask and independently checked with Bareiss.
6. A truncated tie-mask list no longer invalidates completed exhaustive work.
   Bridge endpoints are reconstructed and checked independently.
7. Tie accounting distinguishes engine-mask occurrences from rerooted logical
   masks and excludes logical mask zero explicitly.
8. The runner no longer depends on macOS-specific `/usr/bin/time -lp`.

## Interpretation

The recommendation was executed as intended: use QD outputs to choose
separated high-quality parents, use exact adjugate/pair information to build
small affine neighborhoods, and shard exhaustive evaluation. The method is
fast and auditable, but neither the broad core tranche nor the completely
fresh 32-bit connector moved the frontier.

The next affine tranche should avoid replaying known frontier endpoints.
The most defensible unused candidates are the subfrontier pairs 018↔022
(28-entry bridge) and 018↔024 (22-entry bridge plus transverse repairs).
They should first pass the same exact affine-intersection audit.

## Pinned hashes

```text
1ae128dcd18acc1ccb6051c42bcae59d96afdbf2252ed87f291335fd912bc362  research/qd_cube_campaign.py at plan time
6a39274319aa3596737f60bebcfbb8da22efaf3acc8a2ee269b90ba09a6da9a3  research/fast_cube_partition32.py at connector run time
4632f5a24f2a2e9ee106cfbe260ee0e0bd9424099c3dfb3f6aaa9a1fafe017fb  research/fast_principal_cube.cpp
9e23bc5ce38eba65d0366bff35b89e91a0d6dbc0255566530fda2f05e9279ecd  build/research/fast_principal_cube_lnps
5bb9da58930b525250799b5a0142c1348435cc8cc0174d29abb3dfea3f6d20fb  QD-selected aggregate report
4d5d85fb4ed4e8ef30cb7242015a942adad52aca9502024c985e843526e7b693  aligned-connector aggregate report
12fe6eefafed40dffab5b8f45bd9b66354222e18165ba466058f8721492e0f9d  research/affine_gf2_audit.py
cafb757a90ae0ce4798fb12159808029a7cb83f5b0e8ea0a45c409b4a5ff47e4  aligned-connector affine audit report
5dce5c286ad6233ce92d51eb5d298b51712bc98204bc07879cbbc45aecbe33a1  research/fast_cube_connector_provenance.py
ca3613f55e83ea84aaeaa9d3cc607a822bec9049a4611ebb33dd6b1ae1cb19d3  aligned-connector provenance sidecar
94a7f563f702e1750ed5459c55e6be977764fc22c37bbc07c39bced1355ba925  aligned-connector 294-file artifact inventory
```

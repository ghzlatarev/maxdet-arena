# Exact top-K beam LNPS and two deep 32-cubes

Date: 2026-07-28

## Outcome

No matrix exceeded the verified order-23 frontier

```text
2,779,447,296,000,000
```

The new work evaluated 110 overlapping 27-bit LNPS cubes and two complete,
individually disjoint 32-bit cubes.  Together with the preceding 200-cube
batch, the exact FPM campaign now accounts for:

```text
evaluated top-level affine regions:       312
27-bit FPM leaf runs:                     374
assignment-visits:             50,197,430,272
engine seconds:                       789.571092
summed wall seconds:                  840.488039
strict promotions:                             0
real-run zero-pivot corrections:               0
best strictly subfrontier score: 2,726,756,352,000,000
```

“Assignment-visits” is deliberate.  The 200 initial cubes and 110 LNPS cubes
can overlap, so their states are not claimed to be unique matrices.  Each
individual 32-cube is a proven disjoint partition of exactly `2^32`
assignments, but the two 32-cubes and the earlier cubes are not claimed
mutually disjoint.

## Evaluator extension

`fast_principal_cube.cpp` now has bounded deterministic capture:

- `--top-k N`: best nonzero, strictly subfrontier masks globally;
- `--top-k-per-weight N`: the same separately for each Hamming weight;
- `--top-k-output-dir DIR`: optional Bareiss-verified global top-K artifacts.

Ranking is exact absolute determinant descending, then mask ascending.
The default path with all three options omitted is unchanged.

Validation completed:

```text
strict C++20 build with -Wall -Wextra -Wshadow -Werror -pedantic: pass
ASan + UBSan self-test:                                      pass
256 principal-minor differential cases:                     pass
72 complete entry-flip cube/Bareiss differential cases:     pass
top-K off/on known-24-cube output and score identity:        pass
39 trusted unittest cases:                                  pass
```

The known 24-bit cube took `0.227631 s` with capture disabled and `0.269466 s`
with global/per-weight capture enabled; both returned byte-identical best
matrices and identical search counters.

## Iterative LNPS pilot

`fast_cube_lnps.py` used exact frontier roots from H0, H1, and H2.  Before
search it checked both directions of:

- the 12-entry H0↔H1 A0 switch; and
- the 12-entry QUBO↔Orrick-reference switch.

Generation 0 scanned six 27-bit supports per root.  Later generations used a
six-state beam per lineage and exact 6/9/12/15 support overlaps.  Every child
support omitted an active incoming flip, making its immediate parent
unreachable.  Exact one-flip, pair-score, and pair-synergy ranks supplied
repair entries.  Beam states were sign-normalized for identity, selected with
Hamming diversity, and independently Bareiss-verified.

The 315-second soft cap produced:

```text
generation 0: 18 / 18 cubes complete
generation 1: 54 / 54 cubes complete
generation 2: 38 / 54 cubes complete
total:       110 / 126 cubes complete
unrun:        16 planned/reserved cubes
assignment-visits: 14,763,950,080
retained weight-4..23 candidate records: 4,400
frontier re-entry masks: 38
strict promotions: 0
```

The old aggregate field `affine_fingerprints_total_seen=326` meant 200 prior
fingerprints plus all 126 newly planned/reserved fingerprints.  It did **not**
mean 326 evaluated cubes.  The corrected report says:

```text
affine_fingerprints_evaluated:                  110
affine_fingerprints_planned_or_reserved_total: 326
```

The original report bytes are preserved, and the corrected report records
that the change was metadata-only.

Evidence:

```text
runs/direct-search/fast-principal-cube/
  lnps-beam-h012-20260728/aggregate-report.json
  lnps-beam-h012-20260728/aggregate-report.pre-metadata-correction.json
  lnps-beam-h012-20260728/provenance/fast_cube_lnps.run-source.py
```

## Frontier tie audit

Five sign-normalized frontier orientations were preserved.  All five passed
both exact Bareiss and `./arena verify`.

- Two sign-normalized-match members of the known 24-state neutral network.
- Three were sign-normalized-new relative to that network plus the QUBO seed.

The latter wording is not an H-equivalence claim.  A full audit with pinned
`pynauty==2.8.8.1` classified the three remaining orientations:

- two in the already represented reference H/HT class `db2cddf4...`;
- one in the already represented QUBO seventh class `b64c3309...`.

The combined local corpus remains seven H-classes, seven HT-classes, and one
normalized row-Gram class.  The LNPS pilot added no local H/HT class.

Evidence:

```text
runs/direct-search/fast-principal-cube/
  lnps-beam-h012-20260728/tie-audit/normalized-comparison.json
  lnps-beam-h012-20260728/tie-audit/h-equivalence-audit.json
```

## Deep 32-bit cubes

Each full support was split by five fixed outer bits into 32 disjoint
27-bit leaves.  The global mapping is

```text
inner_global = engine_mask XOR reroot_xor
global = inner_global OR (outer_mask << 27)
```

All 64 leaf bases were nonsingular, so every recorded `reroot_xor` was zero.
The reroot path remains implemented for future cubes.  For each cube:

- 32 unique affine-leaf fingerprints were recorded;
- the leaf union is exactly `32 × 2^27 = 2^32`;
- 256 independent random masks reproduced the same matrix through direct
  global application and leaf/XOR application;
- all global top-32 matrices were rechecked by Bareiss; and
- all frontier assignments were rechecked by Bareiss and `./arena verify`.

| Full cube | Full fingerprint | Assignments | Wall | Frontier assignments | Best subfrontier |
|---|---|---:|---:|---|---:|
| H2 QUBO↔reference 12-bridge + 15 transverse + 5 pair rescue | `d1f3ccc1ca26d31085e1159aa564a93172d609a14b25abd73607a05ad845acb0` | 4,294,967,296 | 67.896792 s | global masks `0`, `4095` | 2,726,756,352,000,000 |
| H0↔H1 A0 12-bridge + 15 transverse + 5 pair rescue | `dddf4fd35c3c80d84d52a155b38c03b382f0e710b655a419a89443b5b1034f98` | 4,294,967,296 | 66.632376 s | global masks `0`, `4095` | 2,726,756,352,000,000 |

The H0 full fingerprint was checked against the initial 200 fingerprints, all
126 LNPS planned/reserved fingerprints, and the H2 full fingerprint.  It was
new to that set.  “New fingerprint” means a new affine cube, not a new matrix
class.

Evidence:

```text
runs/direct-search/fast-principal-cube/
  partition32-h2-bridge-20260728/aggregate-report.json
  partition32-h2-bridge-20260728/provenance.json
  partition32-h0-a0-20260728/aggregate-report.json
  partition32-h0-a0-20260728/manifest.json
```

## Safe aggregate accounting

| Campaign | Evaluated regions/leaves | Assignment-visits | Engine seconds | Wall seconds |
|---|---:|---:|---:|---:|
| Initial mixed H0/H1/H2 batch | 200 27-bit cubes | 26,843,545,600 | 383.359945 | 390.922356 |
| Iterative LNPS | 110 of 126 27-bit cubes | 14,763,950,080 | 273.839666 | 315.036515 |
| H2 deep cube | 32 leaves = one full 32-cube | 4,294,967,296 | 66.845186 | 67.896792 |
| H0 deep cube | 32 leaves = one full 32-cube | 4,294,967,296 | 65.526295 | 66.632376 |
| **Total** | **374 FPM leaf runs** | **50,197,430,272** | **789.571092** | **840.488039** |

These totals are exact visit counts and summed timings.  They are not a
global proof about all order-23 sign matrices and not a count of unique
matrices.

## Provenance hashes

```text
4632f5a24f2a2e9ee106cfbe260ee0e0bd9424099c3dfb3f6aaa9a1fafe017fb  research/fast_principal_cube.cpp
9e23bc5ce38eba65d0366bff35b89e91a0d6dbc0255566530fda2f05e9279ecd  build/research/fast_principal_cube_lnps
2a9274e34ddfa860328ae2c7b5184eade7ea6818c0ddbe9429360d94f5c41f63  research/fast_cube_lnps.py
00bacd869ec59512879222683654e30a29dd04e982034d1ada6a072698bf0eb8  LNPS reconstructed run-source snapshot
410419286ef4b2a92c1508600344c6a7f4efa020e96756453b904ef30aee3de8  research/fast_cube_partition32.py
8d1462d2b702c9cb28828fdcf59ee45fa513e5c2f5c871320484ac98e8488400  H2 partition driver at run time
2ad60931d99d5d5874b7969dcaa46ee05f87e088e95ad1bea6256a52f75e3  initial-200 aggregate
2f0a585894e6cd9507d2d3ac72ff1ddd5d7c9d15483bdcc32372064065f69f2f  corrected LNPS aggregate
8ac48ab4e0217f3942654e7129f0cec2980cde27f1740ac11f7744f5e0655104  original LNPS aggregate
4f746359c520f13383e8f4ad985d952bcf509358e7dc939482a60ed2a4e8b633  H2 32-cube aggregate
235295851c6d5bd3ccadf6e3f04c46c0db45506b4563ca1855911f62c3a6636f  H2 provenance sidecar
7832cbb623c7dd94b0991c5bf31cc9aeedced41baccdcc62ffc1a22e92ae853a  H0 32-cube aggregate
```

The H0 report embeds the driver, FPM source, and FPM binary hashes directly.
The earlier H2 report is unchanged; its sidecar binds the immutable report to
the at-run driver hash and the same FPM source/binary.

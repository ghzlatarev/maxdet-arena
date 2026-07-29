# Exact radius-3 Gram orbit screen

**Run date:** 2026-07-28

## Result

The complete neighborhood obtained from the published 45-edge defect graph
by deleting exactly three present edges and adding exactly three absent edges
contains

```text
C(45,3) C(208,3) = 20,976,452,640
```

labeled graphs. Exact quotienting by the full automorphism group reduces this
to `9,967,496` representatives. Every representative was determinant-screened:

```text
positive determinant representatives       9,967,496
determinants above the squared frontier         1,958
perfect-square determinant representatives         55
perfect squares above the frontier                    0
qualified Hasse/shell route hits                       0
```

The 55 square orbits represent 22,488 labeled graphs and all have roots below
the frontier. Their largest root is `2,759,196,672,000,000`.

The largest Gram determinant in the neighborhood is

```text
7,985,301,658,848,460,800,000,000,000,000
```

at orbit `787708`, of size 12. It is positive definite but not a square. An
independent Python Bareiss replay reproduced it exactly.

Therefore this exact neighborhood contains no Gram candidate for a strict
matrix improvement. This is a local Gram result, not a global MaxDet or
optimality claim.

## Complete orbit argument

The pinned defect graph has one 15-vertex component and two isolated `K4`
components. Its full automorphism group is

```text
((S2 x S2)^3 semidirect S3) x (S4 wr S2)
```

of order

```text
384 x 1,152 = 442,368.
```

`gram_radius3_orbits.cpp` enumerates all 442,368 permutations and validates
that each preserves the exact published edge set. It first takes orbits of
three-edge deletion subsets. For each deletion representative it enumerates
the full stabilizer and takes its orbits on three-edge addition subsets. The
pair-orbit size is the deletion-orbit size times the stabilizer-orbit size.

The resulting 287 deletion orbits produce 9,967,496 pair orbits, whose
multiplicities sum exactly to 20,976,452,640. A separate NetworkX group
enumeration and edge-cycle Burnside implementation independently obtained

```text
fixed-pair sum: 4,409,301,270,528
group order:          442,368
quotient:           9,967,496
remainder:                  0
```

Two materializations, using ordinary sorting and a fixed compare/swap network,
produced byte-identical 303 MiB catalogs.

## Exact screening

Each candidate differs from the fixed Gram by a rank-at-most-12 update.
For each of four proven 31-bit primes, the screen precomputes the fixed Gram
inverse and applies

```text
det(B + P Delta P^T)
  = det(B) det(I + Delta P^T B^-1 P).
```

Centered CRT reconstruction is unique under the explicit Hadamard row-norm
bound. Every radius-1 and radius-2 orbit representative was compared with
full 23-by-23 elimination at all four residues and after CRT reconstruction.
The regressions reproduced:

```text
radius 1: 9,360 labeled cases; 12 squares at 2,743,271,424,000,000
radius 2: 21,312,720 labeled cases; 12 squares at 2,779,447,296,000,000
```

The audit also compared 20,000 deterministic radius-3 samples, checked an
independent fifth prime, full-checked every 4,096th screened orbit, and
full-checked every observed square. Strict-warning and ASan/UBSan self-tests
passed.

Catalog construction took 1,251.840 and 982.446 seconds in the two retained
runs. The final unsharded exact screen took 83.268 seconds. It emits a truthful
empty `gram-radius3-orbits` route snapshot; the Hasse and sign-shell stages
have nothing to analyze because no square root is above the frontier.

## Reproduce

```sh
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wshadow -Werror \
  research/gram_radius3_orbits.cpp \
  -o build/research/gram_radius3_orbits

build/research/gram_radius3_orbits --self-test
build/research/gram_radius3_orbits --count-orbits

build/research/gram_radius3_orbits \
  --build-catalog \
  --catalog runs/direct-search/gram-radius3-orbits-20260728/catalog.tsv \
  --manifest \
    runs/direct-search/gram-radius3-orbits-20260728/catalog-manifest.json

build/research/gram_radius3_orbits \
  --screen \
  --catalog runs/direct-search/gram-radius3-orbits-20260728/catalog.tsv \
  --catalog-sha256 \
    b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed \
  --output runs/direct-search/gram-radius3-orbits-20260728/screen.json \
  --route-snapshot \
    runs/direct-search/gram-radius3-orbits-20260728/route.json
```

Screening can be divided deterministically with `--shard-count N` and
`--shard-index I`; an unsharded run was retained to avoid aggregation
ambiguity. Every screen fails closed unless the caller supplies the exact
catalog SHA-256.

## Frozen hashes

```text
source
d4b0038d46c21b0a7af072a27cd3fdd5b3f29d6cc9ab7cae815d89f3ee60f512

final screening binary
309ebbc66cea402f546b7db6efce454204ab1e48f4bfcfb609b98dd6cc7e051e

catalog (both materializations)
b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed

original / fast catalog manifests
6f0b3e20be711bf0da3eab648eb91ba7caddf73a51f941e307bffaf468921d28
64d4569464a82454a910cee294db7847360e48219a6b08d4ac2676afa2009f25

screen report
401da43687b47f0ecdd40b8db9be44e66791516478c405af1a26b63937698edb

route snapshot
6f5852d685a062c3bb301f1148956b1a6ae8e59f542ed4d0f889d27fc50557f5

provenance sidecar
052b421052e54b576b50f9531c23869e67b8209e6a96139bf42a7ae04e6c46d4
```

The complete machine-readable provenance is
`runs/direct-search/gram-radius3-orbits-20260728/provenance.json`.

# Gram-basin hopper and coronal follow-up

## Outcome

The 2026-07-29 pilot and exploitation waves did **not** improve the current
order-23 determinant frontier:

```text
known frontier       2,779,447,296,000,000
first strict score   2,779,447,300,194,304
strict wins                                  0
```

They did produce three useful results:

1. 220 exact HT-Gram basins outside the eight pilot seed basins were reached.
2. Three high-scoring basins were absent from the frozen 20-representative
   comparison set, and a dedicated arm connected the best of them back to the
   known frontier.
3. Rewriting the normalized Gram determinant separated a promising basin's
   objective into a strong block determinant and a weak coronal factor. That
   gives the next search a more informative two-objective archive.

All determinant claims below use exact integer arithmetic. Gram-basin novelty
is only relative to the explicitly classified local corpus. It is not a
literature-wide novelty, matrix-factor-equivalence, or optimality claim.

## Search machinery

`research/core_gram_basin_hopper.cpp` runs exact moving-parent reactive tabu
search over dephased 22 by 22 binary cores. Its archive is keyed by an
invariant, noncanonical signed-Gram sketch:

- odd-order line-product switching;
- mod-4 normalized Gram labels;
- a global label histogram;
- sorted incident-label histograms; and
- an unordered row/column pair.

A different sketch proves a different signed-Gram orbit. A matching sketch
does not prove equivalence. Exact equivalence claims use the separate pinned
`pynauty==2.8.8.1` classifier in
`research/generalized_gram_basin.py`.

The hardened engine:

- gates every non-seed archive observation with `--gate-all`;
- distinguishes initialization from search-added discoveries;
- checks deadlines and interrupts inside kicks and polishing;
- refuses output aliases and archive collisions;
- writes outputs atomically; and
- checks the final exact determinant before publication.

The CLI receipt path was also hardened after an audit command exposed an
ambiguous `--json` invocation. Receipt outputs must now end in `.json` and
cannot overwrite the input matrix.

## Pilot

Four 120-second arms used 12, 24, 48, and 72-bit coordinated kicks from eight
exactly distinct HT-Gram seed basins.

| Measure | Result |
| --- | ---: |
| Completed epochs | 5,305 |
| Tabu moves | 96,978,140 |
| Exact move-direction evaluations | 51,304,244,521 |
| Exported elites independently verified | 256 / 256 |
| Exact HT-Gram basins among 8 seeds + 256 exports | 228 |
| Export basins absent from the 8 seeds | 220 |
| Strict wins | 0 |

Five export basins absent from the eight seeds exceeded
`2,600,468,480,000,000`. Comparing the best five against a frozen set of 20
high-scoring local representatives left these three locally unseen basins:

| ID | Absolute determinant | Core quotient | Exact HT-Gram key |
| --- | ---: | ---: | --- |
| c1 | 2,638,361,395,200,000 | 629,034,375 | `fcc618bf2ab4…` |
| c4 | 2,610,279,795,916,800 | 622,339,200 | `5b458ededf89…` |
| c5 | 2,602,405,829,017,600 | 620,461,900 | `d742055ae55d…` |

An exact radius-two audit tested the center, all 484 one-bit moves, and all
116,886 two-bit moves around every dephased core. All three are strict local
maxima. Their best two-bit descendants have core quotients 616,335,000,
600,661,600, and 604,804,100 respectively.

## Exploitation wave

The follow-up ran four 180-second, all-gated arms:

| Arm | Seeds | Kick | Best core quotient | Search-added sketches |
| --- | --- | ---: | ---: | ---: |
| c1 | c1 | 12 | 662,671,875 | 1 |
| c4 | c4 | 24 | 638,666,000 | 2 |
| c5 | c5 | 12 | 620,461,900 | 0 |
| mixed | c1, c4, c5, frontier, pair-kick | 48 | 662,671,875 | 3 |

Together the arms completed 9,443 epochs, 141,955,052 tabu moves, and
75,099,217,855 exact move-direction evaluations. Every one of the 14 archive
exports and four final best files passed the independent arena verifier.

The c1-only arm reached a byte-distinct frontier matrix. Exact classification
assigned it to H/HT certificate `de7642266b69…`, already represented in the
larger polar-portal corpus. This is a useful basin-to-frontier transition, not
a new frontier class. The mixed arm also produced one additional exact
HT-Gram basin at `2,601,076,654,080,000` outside the frozen 23-basin
high-score comparison set.

## Matched determinant-only control

While the coronal selector was being implemented, a matched five-minute
control used its planned seed groupings with the original determinant-only
parent policy and a lower all-observation gate of core quotient 600,000,000.

| Measure | Result |
| --- | ---: |
| Completed epochs | 13,296 |
| Tabu moves | 242,809,420 |
| Exact move-direction evaluations | 128,453,216,764 |
| Search-added Gram sketches | 34 |
| Archive exports | 43 |
| Independently verified best + archive files | 47 / 47 |
| Strict wins | 0 |

Exact classification of the 43 exports together with the frozen comparison
set found 23 additional HT-Gram basins. Five exceeded
`2,600,468,480,000,000`. The two strongest were:

| Absolute determinant | Core quotient | Exact HT-Gram key |
| ---: | ---: | --- |
| 2,654,208,000,000,000 | 632,812,500 | `b978ddd68381…` |
| 2,653,274,767,360,000 | 632,590,000 | `3b76bc43eb14…` |

The classifier comparison set contains every previously retained local
HT-Gram representative above the same threshold. These are therefore new
relative to that frozen high-score local corpus, but no literature-wide
novelty is claimed.

## Exact coronal-Pareto wave

The matched `--coronal-pareto` wave used the same five-minute budget, lower
gate, kick sizes, and seed geometry as the control. Parent selection mixed
exact determinant, larger `det(M)`, smaller `kappa`, the current Pareto front,
and exploration. Strict promotion remained determinant-only.

| Measure | Result |
| --- | ---: |
| Completed epochs | 13,199 |
| Tabu moves | 241,138,156 |
| Exact move-direction evaluations | 127,569,066,795 |
| Search-added Gram sketches | 32 |
| Archive exports | 46 |
| Independently verified best + archive files | 50 / 50 |
| Strict wins | 0 |

Relative to every exact basin in the pilot, exploitation, and matched-control
waves, the Pareto wave added 18 exact HT-Gram basins. Seven exceeded the
high-score gate. Its strongest new basin tied the control's best new score,
`2,654,208,000,000,000`, but has a different exact HT-Gram certificate
`c9e8050c93b8…`.

The c4-only Pareto arm improved its matched control result from core quotient
638,666,000 to 640,828,125. The retained exact trade-off chain included:

| Core quotient | `det(M)` | `kappa` |
| ---: | ---: | ---: |
| 638,666,000 | 348,009,506,768,000,000 | 2.827922 |
| 627,004,000 | 319,824,571,136,000,000 | 2.770782 |
| 627,200,000 | 314,703,872,000,000,000 | 2.750000 |
| 640,828,125 | 299,812,213,476,562,500 | 2.630274 |
| 640,828,125 | 291,607,546,289,062,500 | 2.591735 |
| 662,671,875 | 267,647,730,468,750,000 | 2.359284 |

No retained point dominated the frontier on both coronal objectives. The
experiment nevertheless confirmed that the exact decomposition changes
archive behavior and preserves a useful route through realizable matrices.

## A different representation

Let `H` be a normalized order-23 sign matrix, `G = H H^T`, and switch its
lines so every off-diagonal Gram entry is 3 modulo 4. Define

\[
W=\frac{G-(24I-J)}4,\qquad M=6I+W,\qquad
\kappa=\mathbf 1^T M^{-1}\mathbf 1.
\]

Because `G = 4M - J`, the matrix determinant lemma gives

\[
\det(G)=4^{22}\det(M)(4-\kappa).
\]

For the exact core quotient \(q=|\det(H)|/2^{22}\), this becomes

\[
q^2=\det(M)(4-\kappa).
\]

This exposes two objectives that ordinary determinant search multiplies
together too early:

| Basin | `det(M)` | Exact `kappa` | Decimal |
| --- | ---: | ---: | ---: |
| c1 | 236,046,727,089,843,750 | 6,240,741,415 / 2,685,687,206 | 2.323704 |
| c4 | 293,534,443,238,400,000 | 2,962 / 1,105 | 2.680543 |
| c5 | 245,651,411,424,020,000 | 59,763,267,634,447 / 24,565,141,142,402 | 2.432849 |
| frontier | 267,647,730,468,750,000 | 6,455 / 2,736 | 2.359284 |

c4 already has a block determinant 9.672% above the frontier's. Its entire
score deficit is the worse coronal factor. Holding `det(M)` fixed, it crosses
the frontier if `kappa` falls below approximately 2.503978.

The normalized c4 defect graph consists of two structured wheel-like
components. A continuous Gram surgery that replaces one wheel by alternating
triangles points just above the current frontier, but its determinant is not a
square. It therefore cannot be `H H^T` and is only a directional relaxation,
never a candidate matrix.

The implemented `--coronal-pareto` mode stores the sorted row/column pair of
exact `(det(M), kappa)` points. This makes the feature invariant under
transpose even though a row coronal alone is not. Archive survival protects
the exact extremes in core quotient, larger `det(M)`, and smaller `kappa`,
then prefers less-dominated points. Parent draws are 30% determinant quality,
20% larger `det(M)`, 20% smaller `kappa`, 20% from the Pareto front, and 10%
least-used/uniform exploration. Exact determinant promotion is unchanged.

Every seed gets a machine-readable `seed_coronal` JSONL event and a console
line. Exported elites carry both exact coronal points in `manifest.json`.
Self-test fixtures reproduce:

| Fixture | `det(M)` | `kappa` |
| --- | ---: | ---: |
| frontier | 267,647,730,468,750,000 | 6,455 / 2,736 |
| c4 | 293,534,443,238,400,000 | 2,962 / 1,105 |

The coordinates are exact diversity heuristics, not an alternative score or
evidence of a distinct H/HT class.

## Independent next representation

The strongest separate follow-up is exact Hamming-code clique repair. After
the same normalization, rows are 23-bit words and

\[
G_{ij}=23-2d_H(x_i,x_j).
\]

The high-quality `{-1, 3}` Gram stratum is therefore the distance set
`{12, 10}`. Destroying 4 to 10 rows, enumerating all compatible 23-bit words,
and solving the replacement as an exact clique jumps directly between
realizable high-Gram basins. Initial pool probes found only about 20 to 400
compatible words for the useful destroy sizes. A Hadamard-Fischer determinant
bound can safely prune partial cliques.

This is an application of exact clique coding methods and large-neighborhood
destroy/repair, rather than another random bit-kick representation.

The first row and transposed campaigns completed 159 exact destroyed-line
neighborhoods across four seeds. They exhaustively tested 1,400,897,536
oriented 23-bit words and traversed 79,601,172 clique nodes. Exact
Hadamard-Fischer bounds pruned every branch before a full leaf; no seed
improved. Eight deadline-interrupted neighborhoods are explicitly excluded
from the closure count. This is exact local evidence for the recorded subsets,
not an exhaustive radius-eight result.

## Reproduction boundary

Core checks:

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Wshadow -Wpedantic -Werror \
  research/core_gram_basin_hopper.cpp \
  -o build/research/core_gram_basin_hopper
build/research/core_gram_basin_hopper --self-test 10000 --seed 38100
/tmp/maxdet-h-audit/bin/python research/generalized_gram_basin.py --self-test
clang++ -std=c++20 -O3 -I/opt/homebrew/include \
  -Wall -Wextra -Wshadow -Wpedantic -Werror \
  research/hamming_clique_lns.cpp \
  -o build/research/hamming_clique_lns
build/research/hamming_clique_lns --self-test
./arena test
```

A coronal arm adds one flag to an otherwise ordinary basin-hopper run:

```sh
build/research/core_gram_basin_hopper --coronal-pareto \
  --seed-matrix SEED.txt --output best.matrix.txt \
  --archive-dir archive --log run.jsonl --summary summary.json \
  --gate-all --quotient-gate 620000000 --seconds 180
```

Campaign configurations, exact summaries, selected matrices, receipts, and
classification reports are retained under:

```text
runs/gram-basin-hopper-pilot-20260729/
runs/gram-basin-hopper-exploitation-20260729/
runs/gram-basin-hopper-coronal-control-20260729/
runs/gram-basin-hopper-coronal-pareto-20260729/
runs/hamming-clique-lns-wave-20260729/
runs/hamming-clique-lns-transpose-wave-20260729/
```

The public visualization consumes a bounded snapshot generated by:

```sh
python3 tools/update_search_progress.py \
  --run-root runs/gram-basin-hopper-coronal-pareto-20260729 \
  --output public/search-progress.json
```

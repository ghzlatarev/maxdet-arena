# Direct MaxDet search campaign — 2026-07-28

## Status

**No strict improvement has been found.**

The exact arena comparison point is

```text
2,779,447,296,000,000
```

The campaign best is an exact tie at

```text
2,779,447,296,000,000
```

Every order-23 sign determinant is divisible by `2^22 = 4,194,304`.
Consequently the first possible strict score is

```text
2,779,447,300,194,304
```

This is the next necessary score-lattice point, not a matrix found by the
campaign. The current verified frontier remains the lower value above.

Independent H24 searches, exact shell factorization, and a deeper multiflip
beam all produced retained matrices at that score. None is a strict
improvement. A later exact-pair QUBO trust-region pilot produced another tie
which expands the frozen local corpus from six to seven H-equivalence classes;
that is a local-classification result, not a literature-novelty claim. The
strongest retained sub-frontier score is

```text
2,762,200,842,240,000
```

It is `17,246,453,760,000` below the comparison point, a relative gap of
approximately `0.6204994%`. Several searches independently reproduced the
comparison score, but a tie is not an improvement.

All scores called “matrix scores” below were recomputed by the trusted
`./arena verify` path. Values produced by the Gram search are explicitly
labelled **Gram-only square roots**. They are not determinants of known sign
matrices and are not records.

This report makes no world-record, optimality, or Hadamard-inequivalence claim.

## Scope and isolation

This was a direct mathematical search, deliberately isolated from the earlier
Agent Harbor campaign:

- worktree: `/Users/georgizlatarev/maxdet-arena-direct-search`
- branch: `research/direct-search-20260728`
- clean base: `origin/main` at `4206ace9005c`
- original Agent Harbor worktree:
  `/Users/georgizlatarev/maxdet-arena`, branch
  `research/ah-maxdet-campaign-20260726`

No Agent Harbor orchestration or campaign harness was run in this worktree.
Only generally useful lessons, such as requiring exact checkpoints and keeping
search state separate from trusted verification, were retained. The search
engines here are small direct C++ or Python programs under `research/`, with
outputs under the ignored `runs/direct-search/` tree.

Before research edits, the clean base passed:

```sh
./arena test
python3 tools/verify_repository.py
```

The first command reported all 39 dependency-free tests passing, and the
repository verification completed successfully.

## Exact matrix artifacts

| Artifact | Exact score | Raw SHA-256 | Sign-normalized SHA-256 | Receipt SHA-256 |
| --- | ---: | --- | --- | --- |
| Published arena reference, `references/orrick-et-al-2003/matrix.txt` | 2,779,447,296,000,000 | `d134c240811076c9f807b98974ca68fda1e7756d1cf7e7ae72cdc044d8743850` | `c0ea58d361945b20dad78bddb3fd93c0810b762e527e8781e87d6db5e86d993a` | `45578d90f1d0e660ee50aa86d3aecd60c859627fdb1bec264024770ae32a32e5` |
| H24 class-14 late reproduction, `runs/direct-search/best-below/frontier-class14-late-28751.matrix.txt` | 2,779,447,296,000,000 | `a0eb1e054d38cb6328d1a8397f84bb11d6e0345c84703cfc8b820d529673d4d0` | `7ef2379777d0f9515ad463804ba50aba3912810ec7d88e8937f6cdf32c2ee4f5` | `0b7cae6fc5562c1742d7f4ffbc7f10e8d682b7c60b456cc2f3cc63aea4057341` |
| H24 class-9 deletion reproduction, `runs/direct-search/h24-deletion-remaining-classes/elites/class-9-r10-c24.matrix.txt` | 2,779,447,296,000,000 | `3fa9fd308dc5634bc6d8ab02b90b995cfb60f7c2f84f1a26acb19fc5a0e77df1` | `efea293d313742c5e0b72167d59bccb5f25384fd51c87dd65e2634e95a6b40f7` | `6585bb1c53fdf922fe606411451f00a67cfec813c27d6903b961f923d335f62d` |
| Twelve-flip beam tie, `runs/direct-search/multiflip/class9-depth12-tie-replay-29831/tie.matrix.txt` | 2,779,447,296,000,000 | `edbf3e96a40f666f00d4d3bfdd82a8d0316317d94ef134cf9ba5bc3fbf3ff90d` | `6f6d5f23de43f9be7b489ace84847ce455cb214e396a35eaa010acdc25176091` | `8e0d1a3edef44e41a379f2320ea590b411f3aeab4cbeec0141fa25af1cf390e3` |
| Second twelve-flip beam tie, `runs/direct-search/multiflip/depth12-tie-harvest-29851/tie.matrix.txt` | 2,779,447,296,000,000 | `b094145c07472526a430a06e2d54c70652bee236b7846be31483120632ea3a3e` | `cea8ee9d79ec63f39ec901da1387698f2dc35b7c7369143a3ad09887dc1ab6a4` | `7f1b6539b58bf29649fb537b62054534e8316425e0408498fc9cd139c1643c9c` |
| Third beam tie, twelve flips from the second, `runs/direct-search/multiflip/tie2-radius12-harvest-29862/tie.matrix.txt` | 2,779,447,296,000,000 | `640813c0b9fbd14e8a8694972edee25008d62126a69fd50c33e1946a4118a6dc` | `96678d284cf5e824c173f1b4cfbc0df902f4f8473a843113ed78482e00ef3d68` | `22ae90ed370a6f9c14793dba3573ea1e80b25fb16f1cf15a18fab13e28fa04e1` |
| Fourth beam tie, twelve flips from the third, `runs/direct-search/multiflip/tie3-radius12-harvest-29911/tie.matrix.txt` | 2,779,447,296,000,000 | `f84f7096365d43cb8370567cb7883d6555f9c20f4ad4d792bde3fd2a035eb49f` | `9e941845741382785fcf3789005c0f3c60a6243ba6c781fd912e964fe77bd412` | `216957f53e2aaa10ee203d8710955aef1f4eb4af8b21f94deaf4ffb8fddc8ff2` |
| QUBO trust-region tie, `runs/qubo-trust-pilot-20260728-seed31003/best-proposal.matrix.txt` | 2,779,447,296,000,000 | `b386a8714a7bba8e2e88796f136e607ae56402f45695b31bc52f87b57ef0ca52` | `e8cd87169219295f3e4bd5bb3eb3b9b9fa6f17e1144589fa4f46a48fd3b74e84` | `dcf04ab5afe78c11b917a730fc13118b50dcbc23f6fe581ca66bc4ab21c65f7a` |
| Radius-24 sphere tie, `runs/direct-search/hamming-sphere-long/tie1-r24-29901/tie.matrix.txt` | 2,779,447,296,000,000 | `6e5f51f4cc7055b9ea32f8987a26e1c8d70c0fb76c25743f324215b9f9d4a410` | `1dec35ac64b529715dafc13384d302565d4ddc347eba294fc173857e2834863e` | `55e6e6171ef0e270cc5a9de7e346db2439acbb160545eca35b086ebbd1a8bf4a` |
| Radius-32 sphere tie, `runs/direct-search/hamming-sphere-long/class9-r32-29902/tie.matrix.txt` | 2,779,447,296,000,000 | `eb29a4613abfeebf2f8a4922b4cc60bbde4bfc5cf54aaf4567dc863721625691` | `c49b23c6fb6785715110141fd044beff85b083df9f571b2ac1dcc989a979655e` | `3b19d580f34510cfdc6232c0e7d14d3cfabadb69c665edf58576f6fabd45b288` |
| Outside-cycle beam tie from the radius-32 factor, `runs/direct-search/multiflip/sphere32-radius12-harvest-29940/tie.matrix.txt` | 2,779,447,296,000,000 | `85a63a49240a8e09abf4cfe9ce4d3dd405e7ab3206b2aaedfc4c14cf2fe13fbc` | `c4336ccad2558ab4f9c653102877f326955ba680ccf38d3561a72072853b51fd` | `ec4ac8cac974c1cb8ff3f0f15638386b2d044978d3945d0c02c6a6c5375fe894` |
| Exact shell-MILP factor, `runs/direct-search/gram-shell-milp-reference-seed29762.matrix.txt` | 2,779,447,296,000,000 | `00433daf0bbbd2e0413c1c3951d80729fbb5fd71dddf46847549a04732d579e1` | `43a705dcb49faece09af013e7cee98a1a095f06f00b7cfddbeaf16ee001cf737` | `608a9407b0839fbb2d0268bdbcdcf4346754872f8cd30d6d5b8e5416db670fd4` |
| Pre-April-2003 archived factor, `runs/direct-search/reference-data/orrick-pre-april2003.matrix.txt` | 2,779,447,296,000,000 | `b6a0801d0d1459ec262ae4e640c1b94b4baa3c44ea4b20aafbd1f8ae62f82be1` | `09a498bbef2d992ebf1fda856006fe20743993a4bc3f1998588d38790b87c202` | `4017a943b67e9120a4a539c1c711c955710b72f44fd715e0ab7d33c48b6845c1` |
| New class-14 near-frontier matrix, `runs/direct-search/audits/class14-near-frontier-audit3-29507.matrix.txt` | 2,762,200,842,240,000 | `57bb0839fe1831f79449e8578921a50ac458817ad3585634499f49d697592b8f` | `8b92291df23309f01058ec0b640da3fc17529d0867acbf937ebaa82a3ee05263` | `132c2e8861438bffe3481453f2bfc085a9641797efdc3448eddab577271691a8` |
| Alternate class-42 near-frontier factor, `runs/direct-search/cross/class42-alt-near-29543.matrix.txt` | 2,762,200,842,240,000 | `a2fb78b1cf18d9af8eee5766949ff71b0db1796259ec5818c0ea8568551af89b` | `ed5d307bf73f77a7d39427b60cbd26f611fdb9431da749e37ec22ea045d1ea89` | `1eee216680d6da7df49598e9d7d2d2d6f8d838ab3fc9de85b5f545d9a0af1026` |

The sign-normalized hash removes row and column sign choices only. It does not
canonicalize row or column permutations and therefore does not prove that any
two rows in this table are Hadamard-inequivalent.

`research/h_equivalence_audit.py` subsequently performed a full pivoted,
dephased bipartite-nauty audit of the local frontier corpus. At the pinned
input snapshot, 55 paths contained 41 byte-distinct matrices but only six
H-classes, six classes after also allowing transposition, and one normalized
row-Gram class. The class counts were `28, 5, 4, 2, 1, 1`. This is a
classification of the frozen local files only; Orrick reported at least 14
order-23 record classes, so the local corpus is incomplete. Method,
certificates, exact score-lattice arithmetic, and reproduction instructions
are in `research/FRONTIER_STRUCTURE_20260728.md`.

Recheck any retained matrix with:

```sh
./arena verify PATH/TO/matrix.txt
```

## Historical factor recovered

An Internet Archive capture of the order-23 MaxDet page from before April 2003
contains a different exact factor at the same comparison score:

https://web.archive.org/web/20030308034054id_/http://www.indiana.edu:80/~maxdet/d23.html

A later capture contains the factor now retained in the arena repository:

https://web.archive.org/web/20130426172738id_/http://www.indiana.edu/~maxdet/d23.html

The archived factor differs from the repository reference in 118 entries
across 20 row positions and has different raw and sign-normalized hashes. Its
row Gram and column Gram are exactly the same as those of the repository
reference. These observations establish a second exact factor representation,
not full inequivalence. In particular, the identity row mapping failing to be
a signed column permutation is not a complete equivalence test.

Orrick later reported finding 14 distinct order-23 record matrices by gradient
ascent; that historical result provides context but does not turn the simple
hash comparison above into a classification proof:

https://arxiv.org/abs/math/0511141

The later exact H-equivalence audit does classify the two archived factors:
they belong to the same H-class and remain in one class after allowing
transposition. A dedicated recovery audit of the two arXiv source archives,
the distinct Indiana page captures and site index, Brent's current and
archived indexes, and the Mendeley collection recovered no additional
order-23 factor. This does not contradict the published count; the checked
sources preserve the count but not the other matrix bytes or reproducible
generation material. The source-bound result is documented in
`research/PUBLISHED_ORDER23_CLASS_RECOVERY_20260728.md`.

## Hadamard-24 seeding

All 60 order-24 Hadamard representatives in the downloaded Mendeley archive
were parsed. For each representative, the seed generator first checked
`HH^T = 24I` exactly, then deleted one row and one column to produce a
23-by-23 sign matrix. Every direct minor had determinant magnitude `24^11 =
1,521,681,143,169,024`.

Dataset:

https://data.mendeley.com/datasets/hzf94h43c5/1

Representative extraction:

```sh
python3 research/hadamard24_seed.py \
  --input runs/direct-search/reference-data/mendeley-hzf94h43c5-v1/Text/Hadamard_24.txt \
  --class 51 --delete-row 1 --delete-column 1 \
  --output runs/direct-search/h24-seeds/class-51.matrix.txt
```

### Sixty-class reactive screen

Each class received a reproducible 10-second reactive-tabu screen. In total,
the screen completed 60/60 runs, 600.013 search-seconds, 156,028,138 proposed
single-entry moves, and 2,449,191 exact checks.

The leading screen results were:

| H24 class | Seed | Exact best |
| ---: | ---: | ---: |
| 14 | 28614 | 2,638,361,395,200,000 |
| 51 | 28651 | 2,630,667,468,800,000 |
| 1 | 28601 | 2,597,023,618,433,024 |
| 42 | 28642 | 2,596,405,248,000,000 |
| 54 | 28654 | 2,593,262,665,728,000 |
| 13 | 28613 | 2,571,436,032,000,000 |

Logs and checkpoints are under:

```text
runs/direct-search/h24-reactive/class-N.jsonl
runs/direct-search/h24-reactive/class-N.matrix.txt
```

### Long H24-derived searches

Reactive tabu from class 51, seed `28752`, reached the exact comparison score
after 111,684,142 moves and 422.630 seconds. A second class-51-derived run,
seed `28755`, reached the same score after 246,279,276 moves and 912.658
seconds.

The class-14 run, seed `28751`, first reached
`2,762,200,842,240,000` after 801,829,617 moves, then reached the exact
comparison score after 1,492,431,317 moves and 5,452.995 seconds. Its complete
7,200-second run covered:

```text
moves:         1,894,540,346
exact checks:     29,602,229
best score:    2,779,447,296,000,000
```

The adjacent log is
`runs/direct-search/h24-long/class-14-28751.jsonl`. Subsequent searches from
this checkpoint have not produced a strict improvement. The retained late
factor is listed in the exact-artifact table.

An alternate deletion from H24 class 42, seed `29542`, reached the same exact
near-frontier score after 33,970,924 moves and 150.907 seconds. It has
different raw and sign-normalized hashes, and it independently completed the
full pair and cross audits described below. A signed row/column Gram
comparison is consistent with the two near-frontier factors being equivalent;
no full Hadamard-inequivalence claim is made.

### Complete 60-class deletion sweep

Every one of the `60 × 24 × 24 = 34,560` row/column deletion pivots was then
used as a short reactive-tabu seed. The work was staged with different bounded
budgets:

| Classes | Pivots | Seconds per pivot | Moves | Exact checks | Best |
| --- | ---: | ---: | ---: | ---: | ---: |
| 14, 42, 51 | 1,728 | 0.25 | 130,881,009 | 2,355,467 | 2,694,905,856,000,000 |
| 19, 30, 35, 36, 38, 39 | 3,456 | 0.10 | 98,988,517 | 2,119,448 | 2,580,524,128,272,384 |
| Remaining 51 classes | 29,376 | 0.05 | 330,843,073 | 10,178,814 | 2,779,447,296,000,000 |
| **Total** | **34,560** | — | **560,712,599** | **14,653,729** | **2,779,447,296,000,000** |

The sweep completed every scheduled pivot. It found one exact frontier tie,
from class 9 with row 10 and column 24 deleted, and no strict improvement.
Its exact matrix is listed above. Authoritative summaries are:

```text
runs/direct-search/h24-deletion-14-42-51/summary.json
runs/direct-search/h24-deletion-high-diversity/summary.json
runs/direct-search/h24-deletion-remaining-classes/summary.json
```

The 15 strongest distinct class leaders then received 120 seconds each:
523,030,536 additional moves and 8,172,466 exact checks. Their best was
`2,687,827,968,000,000`, with no tie or win. Two further 1,200-second runs
from classes 6 and 10 added 721,428,227 moves and 11,272,336 exact checks.
Both reached `2,762,200,842,240,000`; neither reached the frontier.

This is complete coverage of the 34,560 direct minors in the pinned
60-representative catalog, followed by bounded local searches. It is not an
exhaustive search of every matrix reachable from those seeds, and it makes no
new claim about the catalog's equivalence classification.

## Direct matrix search methods

All floating-point inverse or cofactor work in these engines is proposal-only.
Best-checkpoint promotion and reported determinant comparisons use exact
integer arithmetic; retained outputs are independently passed to
`./arena verify`.

### Reactive single-entry tabu

`research/reactive_tabu.cpp` ranks all 529 single-entry flips, permits
controlled downhill moves through a reactive tabu tenure, and periodically
rebuilds its inverse. Selected stopped runs alone account for 2,760,408,408
moves and 43,131,583 exact checks over 10,740.582 search-seconds:

| Start | Seed | Moves | Exact checks | Best |
| --- | ---: | ---: | ---: | ---: |
| Published reference | 28302 | 735,735,605 | 11,495,950 | 2,779,447,296,000,000 |
| Pre-April-2003 factor | 28907 | 461,070,976 | 7,204,275 | 2,779,447,296,000,000 |
| Alternate class-51 frontier | 28905 | 270,258,088 | 4,222,825 | 2,779,447,296,000,000 |
| H24 class 51 | 28752 | 886,628,462 | 13,853,597 | 2,779,447,296,000,000 |
| H24 class 1 | 29201 | 406,715,277 | 6,354,936 | 2,630,667,468,800,000 |

This subtotal excludes active or later-rotated jobs and the separate 60-class
screen. None produced a strict improvement.

Seven later 1,800-second campaigns covered four shell-MILP factors and the
first three neutral-chain ties:

| Start family | Runs | Moves | Exact checks | Strict gains |
| --- | ---: | ---: | ---: | ---: |
| Shell-MILP frontier factors | 4 | 2,248,127,038 | 35,127,150 | 0 |
| Neutral-chain frontier ties | 3 | 1,738,466,669 | 27,163,661 | 0 |
| **Total** | **7** | **3,986,593,707** | **62,290,811** | **0** |

All seven completed their full wall-clock bounds and retained exact frontier
checkpoints. These runs tested wider stochastic basins around raw-diverse
factors; the later H-equivalence audit shows that raw diversity alone should
not be mistaken for seven independent design classes.

Representative command:

```sh
build/research/reactive_tabu \
  --start runs/direct-search/h24-reactive/class-14.matrix.txt \
  --output runs/direct-search/h24-long/class-14-28751.matrix.txt \
  --log runs/direct-search/h24-long/class-14-28751.jsonl \
  --seed 28751 --seconds 7200 --heartbeat-seconds 60
```

### Exact two-line replacement

`research/pair_search.cpp` jointly optimizes every unordered pair of rows and
then every unordered pair of columns. It fixes one redundant whole-line sign
and enumerates `2^22` assignments per pair. One complete pass is therefore:

```text
506 line pairs × 4,194,304 assignments = 2,122,317,824 assignments
```

Five direct-search passes completed without a strict improvement:

| Starting score | Seed | Assignments | Log |
| ---: | ---: | ---: | --- |
| 2,638,361,395,200,000 | 28714 | 2,122,317,824 | `runs/direct-search/h24-refine/class-14-pair.jsonl` |
| 2,687,827,968,000,000 | 28753 | 2,122,317,824 | `runs/direct-search/h24-refine/class-51-best-pair.jsonl` |
| 2,694,905,856,000,000 | 28756 | 2,122,317,824 | `runs/direct-search/h24-refine/class-51-2694905-pair.jsonl` |
| 2,762,200,842,240,000 | 29502 | 2,122,317,824 | `runs/direct-search/h24-refine/class14-near-frontier-pair-29502.jsonl` |
| 2,762,200,842,240,000 (alternate factor) | 29544 | 2,122,317,824 | `runs/direct-search/h24-refine/class42-alt-near-pair-29544.jsonl` |

Total completed assignments in these five campaign passes:
`10,611,589,120`. These are local statements for the specified two-line
replacement neighborhood only.

Three raw-distinct neutral-chain frontier matrices later received the same
complete audit, adding `6,366,953,472` exact assignments with no gain. The
three retained outputs independently verify at the unchanged frontier.

The engine also supports entry kicks followed by several settling passes.
Keeping a below-frontier working state separate from the monotonic best output
was essential: otherwise a promising kicked basin could be searched
incorrectly or overwrite a stronger checkpoint.

### Exact row-column cross replacement

`research/cross_search.cpp` optimizes the 45 entries in one selected row and
one selected column jointly. Removing redundant signs reduces each cross to
`2^21` assignments. A full 23-by-23 pass is:

```text
529 row/column crosses × 2,097,152 assignments
  = 1,109,393,408 assignments
```

The following full passes found no strict improvement:

| Start | Seed | Assignments | Elapsed | Log |
| --- | ---: | ---: | ---: | --- |
| Published reference | 29101 | 1,109,393,408 | 44.444 s | `runs/direct-search/cross/reference-audit-29101.jsonl` |
| Pre-April-2003 factor | 29103 | 1,109,393,408 | 45.845 s | `runs/direct-search/cross/orrick-pre-april2003-audit-29103.jsonl` |
| Alternate class-51 frontier | 29104 | 1,109,393,408 | 44.240 s | `runs/direct-search/cross/frontier-class51-audit-29104.jsonl` |
| New class-14 near-frontier matrix | 29501 | 1,109,393,408 | 48.091 s | `runs/direct-search/cross/class14-near-frontier-29501.jsonl` |
| Alternate class-42 near-frontier factor | 29543 | 1,109,393,408 | 44.513 s | `runs/direct-search/cross/class42-alt-near-29543.jsonl` |

The three frontier factors are therefore exact cross-local optima under this
coordinate parameterization after `3,328,180,224` assignments. Including the
two near-frontier factors, the five completed passes cover
`5,546,967,040` assignments.

The same three neutral-chain matrices also completed one cross pass each,
adding another `3,328,180,224` exact assignments with no gain.

A harness quirk was found here: equal-score elites were initially discarded,
and kicks restarted from the original matrix. Equal-score elite rotation was
added, but only for matrices not related by row/column sign changes; merely
cycling sign-equivalent representatives adds no search diversity.

### Coordinated triple-line search

`research/triple_line_search.cpp` perturbs two rows or columns by one to four
entries and replaces a third line with the exact determinant-maximizing
cofactor-sign pattern. It alternates orientations, accepts controlled downhill
moves, and periodically kicks from the exact incumbent.

Five completed/stopped exploratory runs covered 222,933,804 proposals,
222,990,882 exact checks, 54,424 restarts, and 8,711.527 search-seconds. Starts
included the repository reference, the alternate class-51 frontier, the
prelaunch basin, and the class-14 H24 basin. The best score in that subtotal
was a frontier tie; there was no strict improvement.

Principal logs:

```text
runs/direct-search/triple/reference-28903.jsonl
runs/direct-search/triple/frontier-class51-28902.jsonl
runs/direct-search/triple/prelaunch-28901.jsonl
runs/direct-search/triple/prelaunch-28904.jsonl
runs/direct-search/triple/class14-28906.jsonl
runs/direct-search/triple/class14-near-frontier-29505.jsonl
```

Unlike the pair and cross passes, this stochastic triple-line work is not an
exhaustive neighborhood proof.

### Frontier-valley and exact rectangular searches

Four further valley-crossing variants started from exact frontier factors and
mixed small multi-entry changes, whole rectangular kicks, exact optimized-line
settling, cooling, and restarts. Together they covered:

```text
proposals:       128,762,782
exact checks:    128,879,812
restarts:             71,300
search-seconds:         5,100
strict gains:               0
```

The principal logs are
`runs/direct-search/frontier-valley/frontier-{wide,medium,block,block-hot}-*.jsonl`.
These are stochastic negative searches, not exhaustive neighborhood results.

`research/block_screen.cpp` separately completed an exact rectangular audit
around a frontier factor. For every three-row choice, and then every
three-column choice, it flipped the same selected set of 2, 3, 4, or 5
opposite-orientation coordinates:

| Rectangle family | Exact candidates |
| --- | ---: |
| Three rows × 2–4 columns | 19,266,709 |
| Three columns × 2–4 rows | 19,266,709 |
| Three rows × 5 columns | 59,592,379 |
| Three columns × 5 rows | 59,592,379 |
| **Total** | **157,718,176** |

No candidate improved the frontier. The four final logs are
`frontier-block3-{rows,columns}-k{4,5}-*.jsonl` under
`runs/direct-search/frontier-valley/`. This exact statement applies only to
those structured rectangular flips.

### Deeper multiflip and exact neutral-switch searches

`research/multiflip_beam.cpp` was extended from eight to 24 flips and taught
to retain a non-baseline frontier tie separately from the best checkpoint.
Three principal depth campaigns covered:

| Flip depths | Ranked/generated states | Exact checks | Strict gains |
| --- | ---: | ---: | ---: |
| 5–8 | 20,313,793,659 | 18,876,149 | 0 |
| 9–16 | 5,638,427,245 | 23,453,786 | 0 |
| 17–24 | 4,688,810,222 | 16,184,272 | 0 |

The depth-12 harvest discovered a repeatable neutral-switch pattern. Six
pairwise-disjoint 12-entry masks were inferred from four verified ties.
`research/neutral_cycle.py` then checked all `2^6 = 64` compound states by
exact Bareiss determinants. Exactly 12 states met the frontier, none exceeded
it, and the 12 ties form one cycle under single-generator moves. The farthest
cycle state is 72 entries from the class-9 baseline.

An independent sphere-derived tie and its 12-entry neighbor were added as two
more generators. `research/neutral_extension.py` checked all 256 selector
states, which have GF(2) rank eight and produce 256 distinct matrices. It
found exactly 14 frontier states—the original 12-cycle plus a two-edge tail
from the class-9 baseline—and no strict gain. No mixed state joining the tail
to a non-baseline cycle state was another tie.

The outside endpoint then produced a second independent 12-entry neighbor.
`research/neutral_completion.py` derived both missing complements and tested
all four algebraic orientation couplings: 256 selector states, 144 distinct
matrices, and exact Bareiss checks throughout. The unique coupling that forms
a full circuit gives another 12-state cycle and no strict gain. Combining the
two results yields one connected neutral graph with:

```text
verified raw matrices: 24
neutral edges:         25
exact 12-cycles:        2
maximum Hamming distance from class 9: 90
strict gains:            0
```

Ten of the second cycle's ties were not present in the earlier graph. An
independent multiflip run predicted one of them byte for byte before the exact
completion was run.

Finally, `research/neutral_union.py` combined both six-generator systems with
their bridge. The 13 masks have GF(2) rank 13, so all 8,192 selector states
are distinct. Exact enumeration found precisely the same 24 ties and 25
edges, no mixed-cycle tie, and no promotion. Thus the displayed two-cycle
graph is complete inside the full inferred 13-generator compound space.

An exact pivoted bipartite-nauty classification of all 24 retained matrices
splits them into two H-classes of 12 matrices each. Allowing transposition
leaves the same two 12-member classes, while all 24 matrices share one
normalized row-Gram graph class. All 24 twelve-entry cycle edges cross between
the two H-classes; only the 32-entry bridge from selector state `0000` to
`1000` stays within an H-class. Both classes were already represented in the
frozen local corpus. The neutral network therefore exposes a structured move
between inequivalent factors of the same Gram geometry, but it adds neither a
new local H-class nor a strict frontier improvement.

Evidence:

```text
runs/direct-search/neutral-cycle/class9-six-generator-29920/report.json
runs/direct-search/neutral-cycle/class9-eight-generator-29941/report.json
runs/direct-search/neutral-cycle/sphere32-six-generator-29943/report.json
runs/direct-search/neutral-cycle/two-cycle-union-29952/report.json
runs/direct-search/neutral-cycle/two-cycle-union-29952/h-equivalence-audit.json
```

The switch signature reduced a larger exact subproblem to a finite labeled
degree-sequence enumeration. `research/fixed_support_enum.cpp` exhausted
6,508,620 masks for each of two fixed 9-by-9 supports around every one of the
12 cycle states:

```text
complete searches:        24
exact determinants: 156,206,880
frontier ties:             24
strict gains:               0
```

Every search found exactly one neutral mask, and all 24 outputs mapped back
into the same known 12-cycle. The modular score is fully exact: determinant
divisibility, the Hadamard bound, and arithmetic modulo `2^32 - 5` uniquely
recover the signed determinant. There is no floating gate in this result.

These matrices are distinct raw artifacts, not 24 inequivalent designs. The
exact audit above establishes precisely two represented H-classes.

### Fixed-radius Hamming-sphere tabu

`research/hamming_sphere_tabu.cpp` maintains an exact Hamming radius and uses
rank-two inverse updates to score every one-out/one-in exchange. Accepted
states are checked exactly; reactive reverse-move tabu, downhill moves, and
deterministic restarts let the walk persist on a sphere rather than greedily
settling.

After strict and sanitizer calibration, four independent 900-second searches
completed:

| Start | Radius | Swap proposals | Exact accepted states | Frontier ties | Strict gains |
| --- | ---: | ---: | ---: | ---: | ---: |
| Depth-12 tie | 24 | 68,960,691,120 | 5,697,235 | 1 | 0 |
| H24 class 9 | 32 | 72,563,332,321 | 4,565,554 | 1 | 0 |
| Shell-MILP factor | 40 | 75,263,468,265 | 3,849,704 | 0 | 0 |
| Published reference | 48 | 76,493,739,361 | 3,314,756 | 0 | 0 |
| **Total** | — | **293,281,231,067** | **17,427,249** | **2** | **0** |

Both retained tie matrices pass `./arena verify`. The radius-24 tie is one
state of the exact 12-cycle above. The radius-32 tie is outside that cycle and
provided the two-generator tail in the 256-state extension. This is a
stochastic negative search over fixed spheres, not an exhaustive Hamming-ball
audit.

The later exact H-equivalence audit identified two singleton local classes
that lacked fixed-sphere coverage. One 450-second radius-24/32/40 run from
each class added:

```text
swap proposals: 96,633,015,505
exact states:         6,096,543
frontier ties:                2
strict gains:                 0
```

Both raw ties verify exactly and canonicalize back to their respective source
H-class; they add no new local equivalence class. Evidence is
`runs/direct-search/hclass-diversity/audit-and-campaign.json`.

### Exact radius-three Hamming audit

The new class-14 near-frontier matrix also completed an exhaustive check of
every matrix obtained by flipping one, two, or three of its 529 entries:

```text
C(529,1) + C(529,2) + C(529,3)
  = 24,673,089 exact candidates
elapsed: 698.228 seconds
strict improvements: 0
```

The final event was `audit_finished`; this was not a time-limited prefix.
The output independently verified at the unchanged score and hashes listed in
the exact-artifact table.

Evidence:

```text
runs/direct-search/audits/class14-near-frontier-audit3-29507.jsonl
runs/direct-search/audits/class14-near-frontier-audit3-29507.matrix.txt
```

This proves only radius-three Hamming-local optimality of this particular
near-frontier factor.

### Radius-four floating screens

`research/radius4_screen.cpp` visited every four-entry perturbation assigned
to each completed shard. Four independent four-shard campaigns covered the
repository reference, the recovered historical factor, the alternate
class-51 frontier factor, and the class-14 frontier factor:

```text
per start:  C(529,4) = 3,226,076,876 combinations
four starts:           12,904,307,504 combinations
strict improvements:   0
```

The rank-four determinant lemma was evaluated in binary64. Every prediction
within `1,000,000,000,000` of the start score was promoted to exact Bareiss,
and a retained top pool was also checked exactly. The observed calibration
errors were at most one determinant unit.

This is a **complete floating screen**, not an exact radius-four audit. The
wide gate and calibration make a numerical miss unlikely, but there is no
formal rounding-error certificate.

The original snapshots predated input-hash recording. A retrospective
manifest now binds the intended inputs, complete snapshots, append-only logs,
best and research outputs, exact scores, and path-relink evidence:

```text
research/DIRECT_SEARCH_PROVENANCE_20260728.json
SHA-256 1836b477a31be38de1310651a8457d68cad14aceaf23065f031d00642837445d
```

It validates with:

```sh
python3 research/build_direct_search_manifest.py \
  --check research/DIRECT_SEARCH_PROVENANCE_20260728.json
```

The binding is necessarily retrospective. The old radius snapshots did not
record start hashes, and old path logs recorded elite paths and exact scores
but not their hashes. The retained zero-promotion radius outputs are
byte-identical to the intended starts, every current determinant agrees, and
each final snapshot equals its final JSONL record. This strongly reconciles
the retained files, but it cannot prove that an input path was unchanged
between historical execution and manifest creation.

### Signed alignment and path relinking

`research/align_elite.py` reduced the frontier-to-near endpoint distances to
77 and 80 entries under signed row and column permutations, with exact
determinant preservation and replayable transformations. Four principal
`research/path_relink.cpp` campaigns then completed 144 directed endpoint
paths, generated 266,871,421 extensions, and made 26,964,696 exact Bareiss
checks over 1,453.028 search-seconds. No campaign produced a strict
improvement.

After the shell MILP produced additional exact factors, three six-path
campaigns added:

| Endpoint family | Directed paths | Generated extensions | Exact checks |
| --- | ---: | ---: | ---: |
| First MILP-factor star | 6 | 20,012,124 | 2,011,394 |
| Cross-family pairings | 6 | 60,154,911 | 3,598,063 |
| Third-factor star | 6 | 47,375,735 | 3,107,349 |
| **Added total** | **18** | **127,542,770** | **8,716,806** |

Across both tranches, path relinking therefore covered 162 directed paths,
394,414,191 generated extensions, and 35,681,502 exact checks without a
strict improvement.

Path relinking retains an exact-score beam after approximate ranking. It is a
pruned heuristic, not an exhaustive statement about every ordering through
the endpoint differences.

### Exact recombinant cubes and restricted crossover tabu

`research/crossover_tabu.cpp` restricts changes to entries where two exact
endpoints differ. Exact mode Gray-code enumerates every recombinant for a
small difference set. The principal completed cubes were:

| Endpoints | Difference dimension | Exact assignments | Frontier ties | Strict gains |
| --- | ---: | ---: | ---: | ---: |
| Class 9 ↔ first depth-12 tie | 12 | 4,096 | 2 | 0 |
| First ↔ second disjoint depth-12 ties | 24 | 16,777,216 | 3 | 0 |
| Second ↔ third tie | 12 | 4,096 | 2 | 0 |
| Third ↔ fourth tie | 12 | 4,096 | 2 | 0 |

The 24-dimensional cube's third tie is exactly the known class-9 baseline;
there is no fourth neutral combination. Further 12-dimensional edge cubes in
the second neutral circuit likewise contained only their endpoints.

Tabu mode ran for 1,800 seconds between aligned class-9 and shell-MILP
factors whose difference dimension is 130:

```text
moves:        860,764,425
exact checks:  29,776,951
restarts:         105,073
strict gains:           0
```

This is a complete result only for the explicitly enumerated small cubes.
The dimension-130 crossover is a bounded stochastic valley search.

### Exact whole-line hybrids

`research/line_hybrid_screen.cpp` tests every matrix obtained by choosing each
differing whole row, or each differing whole column, from one of two aligned
frontier endpoints. Each orientation is a complete `2^k` enumeration for its
`k` differing lines, and every assignment is checked with exact Bareiss
elimination.

Four distinct endpoint campaigns completed:

| Endpoints | Exact row assignments | Exact column assignments |
| --- | ---: | ---: |
| Reference ↔ aligned class 9 | 4,194,304 | 2,097,152 |
| Reference ↔ aligned MILP factor 29763 | 4,194,304 | 8,388,608 |
| Reference ↔ aligned MILP factor 29764 | 4,194,304 | 8,388,608 |
| Aligned class 9 ↔ aligned factor 29763 | 8,388,608 | 8,388,608 |
| **Distinct-search total** | **20,971,520** | **27,262,976** |

That is 48,234,496 exact hybrid matrices with no promotion. Two campaigns
were rerun after adding a fresh-safe tie archive, bringing executed checks to
77,594,624. The reruns recorded four tie assignments each but only two unique
matrices: exactly the two endpoints. There were no intermediate frontier
ties.

The tie archive requires a nonexistent output directory, writes atomically,
deduplicates by full row-major sign bits, and records the orientation,
assignment, and subset mask. Best, research, and archived-tie checkpoints all
independently pass `./arena verify`.

## Restricted exact Gram screens

For a normalized order-23 sign matrix, the row Gram can be written

```text
G = 24I - J + 4A
```

where `A` is the adjacency matrix of a simple defect graph. A positive
definite Gram with square determinant and the correct divisibility is still
only a necessary condition: an exact `R in {−1,+1}^{23×23}` with `RR^T = G`
must be constructed before its square root becomes a matrix score.

`research/gram_edge_search.py` exhaustively screened three restricted graph
neighborhoods:

| Neighborhood | Exact cases | Perfect-square determinants | Result |
| --- | ---: | ---: | --- |
| One-edge-for-one-edge swaps from the published Gram | 9,360 | 12 | All roots `2,743,271,424,000,000`, below the frontier |
| Add up to four edges to the ideal `K3 + 5K4` graph | 96,741,645 | 0 | No square candidate |
| Two-edge-for-two-edge swaps from the published Gram | 21,312,720 | 12 | All roots exactly `2,779,447,296,000,000`; none above |

The four shards for each large screen completed. Artifacts:

```text
runs/direct-search/reference-swap.json
runs/direct-search/gram-screen-max4-shard0.json
runs/direct-search/gram-screen-max4-shard1.json
runs/direct-search/gram-screen-max4-shard2.json
runs/direct-search/gram-screen-max4-shard3.json
runs/direct-search/gram-double-swap-shard0.json
runs/direct-search/gram-double-swap-shard1.json
runs/direct-search/gram-double-swap-shard2.json
runs/direct-search/gram-double-swap-shard3.json
```

Example:

```sh
python3 research/gram_edge_search.py \
  --mode reference-double-swap \
  --shard-count 4 --shard-index 0 \
  --output runs/direct-search/gram-double-swap-shard0.json
```

### Complete three-connector 12-edge shell

The published Gram lies 12 edges beyond the nonsquare ideal
`K3 disjoint-union 5K4` graph. Its added edges are three `K2,2` connectors.
`research/gram_connector12.cpp` exhaustively enumerates the family in which
the three connectors use three distinct pairs of base blocks.

The labeled family has 12,072,240 configurations. Exact canonicalization
under `S3 × (S4 wr S5)` reduces it to 73 graph orbits. The orbit computation
reproduced the published Gram determinant exactly, then evaluated every orbit
with exact integer arithmetic:

```text
orbits:                         73
positive determinants:         73
determinants above frontier²:    4
perfect-square determinants:     1
above-frontier squares:          0
```

The sole square is the published frontier Gram. All four larger Gram
determinants are nonsquares, so none can be a sign-matrix determinant square.
The largest Gram determinant in this family is
`7,862,828,089,422,643,200,000,000,000,000`, but it is nonsquare. No
Hasse, shell, or factor candidate survives.

The complete report is
`runs/direct-search/gram-connector12-complete-20260728.json`; the exact orbit
argument, group action, and validation are documented in
`research/GRAM_CONNECTOR12_20260728.md`.

`research/gram_connector12_reuse.cpp` then allowed connector decompositions
to reuse a base-block pair while requiring all 12 added edges to remain
distinct. Deduplication and the same exact group action expand the complete
union to 12,838,320 labeled graphs and 113 orbits:

```text
positive determinants:        113
determinants above frontier²:   29
perfect-square determinants:     1
above-frontier squares:           0
```

Again, the only square is the published frontier Gram. The new maximum is
`8,296,580,272,619,520,000,000,000,000,000`, also nonsquare. Evidence is
`research/GRAM_CONNECTOR12_REUSE_20260728.md` and
`runs/direct-search/gram-connector12-reuse-complete-20260728.json`. This
closes the stated union of three edge-disjoint `K2,2` connector
decompositions, not all possible 12-edge augmentations.

### Exploratory near-frontier Gram level screen

After row-sign normalization, the class-14 near-frontier row Gram has
off-diagonal counts

```text
-1: 228
 3:  24
-5:   1
```

An exact exploratory calculation tested all 506 one-entry level changes and
127,512 valid two-entry level changes around this Gram. It found no
above-frontier perfect-square determinant. This calculation was performed in
exact arithmetic, but no durable scripted artifact was retained, so it should
be treated as an exploratory result to reproduce—not as one of the documented
exhaustive proof artifacts above.

## Stochastic Gram search: candidates are Gram-only

`research/gram_tabu.cpp` searches arbitrary normalized defect graphs using
hill, tabu, or annealing moves. A log-determinant rank-two update ranks
neighbors, but promotion gates use an exact CRT determinant, exact square
testing, `2^22` root divisibility, and exact positive-definiteness checks.

The retained seed-29302 snapshot is:

```text
runs/direct-search/gram-tabu-long-29302.json
SHA-256: 44e22eee3b74650fcad599bb0150ee46b5a40e3111a2e033cb60d794d1c2ed2e
```

At that checkpoint it contained 256 stored qualified Gram hits in ten
square-root value classes:

| Gram-only square root | Stored hits | Edge counts represented |
| ---: | ---: | --- |
| 2,887,188,480,000,000 | 10 | 39 |
| 2,891,776,000,000,000 | 20 | 43 |
| 2,898,001,920,000,000 | 1 | 37 |
| 2,900,951,040,000,000 | 2 | 41 |
| 2,909,798,400,000,000 | 18 | 37 |
| 2,913,075,200,000,000 | 12 | 41 |
| 2,929,459,200,000,000 | 4 | 35 |
| 2,936,012,800,000,000 | 8 | 32, 37 |
| 2,946,170,880,000,000 | 162 | 32 |
| 2,949,120,000,000,000 | 19 | 34, 36, 40, 41, 42 |

The largest value is about `6.10455%` above the arena comparison point, but it
is **not a matrix determinant**. It is only the square root of a candidate
Gram determinant.

The first `2,946,170,880,000,000` candidate was independently checked for its
exact determinant, square root, `2^22` divisibility, and positive leading
principal minors. Its root factors as
`2^22 * 3^5 * 5^7 * 37`. These checks qualified it for decomposition work;
they did not construct a sign factor.

Representative search:

```sh
build/research/gram_tabu \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/direct-search/gram-tabu-long-29302.json \
  --mode tabu --seed 29302 --seconds 3600 \
  --kick-size 8 --tabu-tenure 13 --restart-iterations 5000 \
  --max-stored-hits 256 \
  --heartbeat-seconds 60 --checkpoint-seconds 30
```

Other exploratory snapshots, including seeds `29303` (tabu) and `29304`
(annealing), found additional square-root and graph classes. They remain
unfactored necessary-condition candidates and are not included in the exact
matrix table.

This distinction is historically important: the decomposition literature
reports that order-23 putative Gram matrices above the known value could be
generated while their sign decompositions could not be found.

https://maths-people.anu.edu.au/~brent/pd/maxdet_cbr_10.pdf

## Exact Gram rejection gates

### Rational Hasse obstruction

`research/gram_hasse.py` diagonalizes a candidate over the rationals and
checks its finite Hasse invariants exactly. If an invariant disagrees with the
identity form, no rational matrix `R` can satisfy `G = RR^T`, hence no sign
factor can exist.

On a complete 256-hit seed-29302 snapshot bound to SHA-256
`b8753020b0a60dc43eac8fcb4947952be3b6303877a6270d4620386ccd791ce5`,
the exact result was:

```text
206 rejected by a rational Hasse obstruction
 50 had no Hasse obstruction
  0 errors
```

The 50 no-obstruction hits represented three square-root values:
`2,891,776,000,000,000`, `2,909,798,400,000,000`, and
`2,949,120,000,000,000`. “No obstruction” is not a factor and is not evidence
that a sign factor exists.

The initial `2,946,170,880,000,000` class was rejected at primes 5 and 37.

Evidence:

```text
runs/direct-search/gram-hasse-all-validation.json
runs/direct-search/gram-hasse-validation.json
```

Command:

```sh
python3 research/gram_hasse.py \
  --snapshot runs/direct-search/gram-tabu-long-29302.json \
  --all-hits \
  --output runs/direct-search/gram-hasse-all-validation.json
```

The Hasse report is bound to the exact snapshot bytes named inside that report.
The active Gram search later changed its checkpoint bytes, so indices must not
be transferred between snapshot hashes.

### Exact sign-column shell obstruction

`research/gram_shell_filter.cpp` applies a stronger sign-specific necessary
condition. If `G = RR^T`, every column sign vector `s` of `R` lies on the
exact shell

```text
s^T G^-1 s = 1,
```

and the outer products of 23 such columns sum to `G`. The filter exhausts all
`2^22` sign vectors up to global negation and checks over `F_3` whether `G`
lies in the span of the shell outer products.

First, the known published Gram passed the gate:

```text
shell size:      1,382
span rank:         154
augmented rank:    154
result:          no shell-span obstruction
```

Evidence: `runs/direct-search/gram-shell-reference.json`.

Then all 256 hits in the current seed-29302 snapshot, bound to SHA-256
`44e22eee3b74650fcad599bb0150ee46b5a40e3111a2e033cb60d794d1c2ed2e`,
were rejected exactly:

```text
 43: empty sign-column shell
213: G outside the mod-3 span of shell outer products
256: total rejected
  0: passed
```

Every rejection includes a sparse mod-3 dual certificate. The shell result is
therefore a proof that none of those 256 candidate Grams has a 23-by-23 sign
factor. It is not merely a failed heuristic decomposition attempt.

Evidence:

```text
runs/direct-search/gram-shell-all-256.json
snapshot SHA-256:
  44e22eee3b74650fcad599bb0150ee46b5a40e3111a2e033cb60d794d1c2ed2e
report SHA-256:
  a3942fdb266ff131089918f622af8b4f51852442e22b2c7fec217d9b12ae685c
```

Command:

```sh
build/research/gram_shell_filter \
  --snapshot runs/direct-search/gram-tabu-long-29302.json \
  --all-hits \
  --output runs/direct-search/gram-shell-all-256.json
```

The practical lesson is to run Hasse and shell gates before expensive framed
decomposition.

### Exact shell-to-factor MILP

For a Gram that survives the shell obstruction, the shell is more than a
filter: it gives a finite exact factorization problem. The shell exporter now
optionally records every normalized sign mask. `research/gram_shell_milp.py`
then solves for nonnegative integer multiplicities `x_s` satisfying

```text
sum_s x_s s s^T = G.
```

The diagonal equations force 23 selected columns. Off-diagonal equations
force their exact outer-product sum to equal the target Gram. SciPy HiGHS is
used to propose an integer solution, but promotion does not trust the solver:
the script reconstructs the sign matrix, checks `RR^T = G` entry by entry,
recomputes its Bareiss determinant, and only then writes the factor.

Applied to the published Gram's 1,382-vector shell, this produced four exact
frontier factors:

| MILP artifact | Raw matrix SHA-256 | Signed-column-canonical SHA-256 |
| --- | --- | --- |
| `gram-shell-milp-reference-29761.matrix.txt` | `a3ef7c0836fae230f7476ed7f97f15142a08a3b1760124310e82f7b715bff1f1` | `877892e4e96cdd954463bb4b06b30082c7a3ef4819b6a222a55f38d1ff2e1705` |
| `gram-shell-milp-reference-seed29762.matrix.txt` | `00433daf0bbbd2e0413c1c3951d80729fbb5fd71dddf46847549a04732d579e1` | `b2accf433633da6fa9d606e284be02ff6282dc81fb8328f08d73fc8cb73289ae` |
| `gram-shell-milp-reference-seed29768.matrix.txt` | `d7599504a3b21239f62c9b6c9678f2b33e2b676f5ea437e778550954f270a1ea` | `420410b12bd28df26ad5f21ffaa23837e1dead54d73797d04554c8c1746a214a` |
| `gram-shell-milp-reference-seed29773.matrix.txt` | `d161a22e1e13e45c4cddb4794b9c8b2140bac9a6a611a945139f666dd27e0e85` | `7a8826705195b34f32b18b840ccafcd09711aa5175f27632dbc5e6d17fc37a91` |

All four independently pass `./arena verify` at
`2,779,447,296,000,000`. Their signed-column-canonical hashes are distinct,
which rules out equality under signed column permutation alone. It does not
prove full Hadamard inequivalence. Random objective seeds select alternate
feasible factors; they do not alter the exact constraints.

An empty-shell parser crash found during audit was fixed. Strict-warning,
ASan/UBSan, exact shell-vector, fresh-output, timeout-status, and positive
factor regressions all pass. A failed or timed-out MILP solve is still not an
impossibility proof; only an exactly reconstructed factor is promoted.

## Multilevel Gram campaigns

`research/gram_multilevel.cpp` extends the binary defect-graph search to
off-diagonal levels `{-5,-1,3}` and `{-5,-1,3,7}`. Floating-point
log-determinants guide moves only. Every promoted hit has an exact CRT
determinant above the comparison point, an exact square root divisible by
`2^22`, and an exact positive-definiteness check. A hit is still only a
candidate Gram until a sign factor is constructed.

The completed campaigns below are independent stochastic searches, not
exhaustive classifications. “Unique” means deduplicated within that campaign;
screen totals count exact neighbor evaluations and can revisit a state.

| Seed and constrained levels | Exact screens | Square observations | Stored unique hits | Hasse rejected / no obstruction | Shell result after Hasse |
| --- | ---: | ---: | ---: | ---: | --- |
| `29611`, tabu, `-5=1..2`, `+7=0` | 20,543,787 | 546 | 538 | 381 / 157 | 157 rejected; 0 pass |
| `29612`, anneal, `-5=1..8`, `+7=0` | 20,585,625 | 179 | 79 | 4 / 75 | 75 rejected; 0 pass |
| `29621`, tabu relocation, `-5=1`, `+7=0` | 14,581,848 | 1,170 | 1,080 | 649 / 431 | 431 rejected; 0 pass |
| `29631`, tabu relocation, `-5=2`, `+7=0` | 20,248,580 | 958 | 843 | 611 / 232 | 232 rejected; 0 pass |
| `29641`, tabu relocation, `-5=0..1`, `+7=1` | 28,385,374 | 1,683 | 1,619 | 1,029 / 590 | 590 rejected; 0 pass |
| `29643`, tabu relocation, `-5=2..8`, `+7=1` | 8,891,315 | 0 | 0 | — | — |
| `29642`, tabu relocation, `-5=0..1`, `+7=2` | 15,977,822 | 0 | 0 | — | — |
| `29644`, anneal, `-5=0..1`, `+7=2` | 2,646,623 | 0 | 0 | — | — |
| `29645`, cool anneal, `-5=0..1`, `+7=2` | 5,074,474 | 0 | 0 | — | — |
| **Recorded-row sum** | **136,935,448** | **4,536** | **4,159** | **2,674 / 1,485** | **1,485 rejected; 0 pass** |

The fixed-one-`-5` campaign initially reached its storage cap. A deterministic
replay with a larger cap reproduced the exact stream, retained all 1,080
unique hits, and reported zero omitted unique observations; only that replay
is counted above. Every other row also reported zero omissions. Across the
1,485 Hasse survivors, the exact sign-shell gate found 83 empty shells and
1,402 mod-3 span obstructions. Thus all 4,159 stored hits in the first five
rows have exact individual impossibility certificates. The four later rows
had no hits to gate. This does not turn the stochastic searches themselves
into exhaustive proofs about unvisited Grams.

The later rows used explicit adaptive stopping. Seed `29643` was stopped
cleanly after 333 seconds because 8,891,315 screens had produced zero
above-frontier neighbors. The fixed-two-`+7` tabu run was stopped after 612
seconds to diversify its stagnant basin: it had 135,856 above-frontier
neighbors but no exact squares. A high-temperature annealing pilot added 2,801
above-frontier neighbors before its temperature was calibrated downward. The
final cool anneal added 747,778. Combined, the three main fixed-two-`+7`
trajectories screened 23,698,919 neighbors, of which 886,435 were above the
frontier, without one exact square. A separate 15-second temperature-tuning
smoke is intentionally excluded from the table. These are bounded negative
search results, not equal-runtime comparisons or exhaustive impossibility
proofs.

The largest Gram-only square root in the four-level campaign was
`2,854,748,160,000,000`. All 45 stored Grams in that value class, and all other
stored hits, were rejected before factorization. The four-level shell parser
was separately calibrated on a known sign factor with one `+7` edge; it passed
with shell size 23,595 and rank 227 of 227. A second known factor with two
`+7` edges also passed, with shell size 14,132 and rank 239 of 239.

Authoritative snapshots and SHA-256 bindings:

```text
gram-multilevel-tabu-29611.json
  563c866d59012d0fe42ad4f05edcee8eea16071a01dc8d79c460a4363f0091c4
gram-multilevel-anneal-29612.json
  db65e9b92bbe93f105f1b5c7c85dbb8d14bcc9d80c7ad4ba2a99479819dbde54
gram-multilevel-relocate-replay-29621.json
  bbac5232d44caddc54e976530e8bf6f58767ad8f05cf932e76798a3d14282a09
gram-multilevel-fixed2-29631.json
  797a148f86484f5ac484b6ba17bb1d72e3224aca49dab566db50ed2062774871
gram-multilevel-plus7-29641.json
  508dcbd2f1f9f7e1b8458053285ccad2e68bf8a9437aeecc5816a7afe2b78845
gram-multilevel-plus7-wide-29643.json
  650ecfd98dfc09c7c8f79aebe57feb9fc0bf2b1535eb0628066d65302a58016a
gram-multilevel-plus7x2-29642.json
  e3a03e685c876124e185c7b23797759db3ee698cc2deca4424ff32463a07b07b
gram-multilevel-plus7x2-anneal-29644.json
  9ad2c549b6ce96d6a41433a7977c79e0d29c05040d07ec082d5adb651be94989
gram-multilevel-plus7x2-anneal-cool-29645.json
  506a9755af486b6d71d19a95202d1e140c15e4d1e537ad37e067c3a65fa8d010
```

For rows with stored hits, Hasse and shell reports use the same basename with
`-hasse.json` and `-shell.json`. The fixed-one-`+7` reports have SHA-256
`6906a5b1177f1ac5a952d32cde6cbdf33350cf4b8b19ec37e8519f9e344666c6`
and
`609e44193dc1e1f65ec15d1f6c6e782f5232649a484acc051c799e31f610e17c`,
respectively. The known-factor calibration is
`runs/direct-search/gram-shell-known-plus7-validation-29641-current.json`;
the fixed-two-`+7` calibration is
`runs/direct-search/gram-shell-known-plus7x2-validation-29642.json`.

## Framed Gram decomposer validation

`research/gram_decompose.cpp` implements the framed row-decomposition method:
columns sharing the same current prefix are grouped into frames, and the next
row is represented by bounded integer counts within those frames. Every
emitted leaf is checked exactly for `RR^T = G` and for the target determinant.

### Known-factor regression

Anchored regression runs against the known comparison Gram succeeded:

| Seed | Anchor | Exact factor checks | Factors written | Elapsed |
| ---: | --- | ---: | ---: | ---: |
| 23012 | Published factor | 11 | 10 | 0.091 s |
| 23013 | Historical factor | 11 | 10 | 0.079 s |

All 20 emitted matrices independently passed `./arena verify`.

Logs:

```text
runs/direct-search/gram-factors/seed-23012-reviewed.jsonl
runs/direct-search/gram-factors/seed-23013-reviewed-historical.jsonl
```

These runs validate parsing, frame construction, leaf reconstruction, exact
Gram checking, exact determinant checking, hashing, and atomic output. They do
not measure pure randomized discovery because an anchored known path is
deliberately available.

The emitted column-canonical hash removes signed column permutations only. It
does not establish full Hadamard inequivalence among emitted factors.

### Pure randomized calibration

Two `--no-anchor` calibrations were run for 120 seconds each against the known
factorable published Gram:

| Seed | Restarts completed | Nodes | Bounded assignments | Maximum depth | Factors |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 23121 | 2,138 | 108,083 | 477,153,160 | 18 | 0 |
| 23122 | 2,172 | 108,234 | 480,573,260 | 17 | 0 |

Combined, pure randomized calibration explored 216,317 nodes and 957,726,420
bounded assignments without recovering even the known factor. This is an
important negative calibration: failure to find a factor for a new Gram over a
similar budget is weak evidence and must not be called non-decomposability.
Only an exact obstruction such as the Hasse or shell certificates supports a
rejection.

Logs:

```text
runs/direct-search/gram-decompose-calibration/seed-23121.jsonl
runs/direct-search/gram-decompose-calibration/seed-23122.jsonl
```

Representative pure calibration:

```sh
build/research/gram_decompose \
  --input references/orrick-et-al-2003/matrix.txt \
  --output-dir runs/direct-search/gram-decompose-calibration/seed-23121 \
  --log runs/direct-search/gram-decompose-calibration/seed-23121.jsonl \
  --seed 23121 --seconds 120 --max-solutions 1 \
  --no-anchor --branch-probability 0.3
```

Two earlier 60-second Z3 factor models returned `unknown` at their timeouts.
They are not unsatisfiability results and are not used as evidence here. A
malformed experimental PB invocation is likewise excluded entirely.

Primary descriptions of the framed and randomized decomposition methods:

https://arxiv.org/abs/1112.4160

https://arxiv.org/abs/1112.4671

## Audit interpretation

The completed pair, cross, and radius-three passes prove only local statements
in their precisely defined neighborhoods. The radius-four result is a complete
floating screen without a rounding certificate. The stochastic reactive,
triple-line, path-relink, crossover-tabu, and fixed-sphere runs provide search
coverage but no proof of local or global optimality. The neutral 13-generator
union, its two fixed supports, and the stated three-connector Gram family are
exactly complete only in their explicitly defined finite spaces. Gram square,
divisibility, and positive-definiteness gates are necessary conditions only.
Hasse rejection rules out rational decomposition; the sign-column shell
rejection rules out the required sign decomposition.

The following are specifically **not** established:

- that the arena comparison point is globally optimal;
- that it is the current world record outside the pinned arena context;
- that the 24 raw neutral-network matrices are 24 Hadamard-inequivalent
  designs—the exact local audit instead finds two H-classes of 12;
- that the six-class local H audit covers all published frontier classes;
- that a Gram-only square root is attainable by any sign matrix;
- that a timed-out or unsuccessful randomized decomposition is impossible.

## Current conclusion

The campaign has not beaten `2,779,447,296,000,000`.
The value `2,779,447,300,194,304` is the first possible score-lattice point
above it, not a discovered matrix.

It did produce:

1. complete direct-minor coverage of all 60 pinned H24 classes, followed by
   bounded refinements and several exact frontier reproductions;
2. two connected 12-state neutral cycles: 24 verified raw frontier matrices
   split into two H-classes of 12, 25 neutral edges, and a maximum distance of
   90 from the class-9 baseline;
3. exact closure of all 8,192 states in the inferred 13-generator union and
   156,206,880 fixed-support degree masks, with no strict gain;
4. a full local H-equivalence audit reducing 41 byte-distinct frontier
   matrices to six H-classes and one normalized Gram class;
5. complete exact pair, cross, recombinant-cube, rectangular, line-hybrid,
   and stated Gram-neighborhood negatives, plus much larger bounded
   stochastic searches;
6. a complete 113-orbit audit of the edge-disjoint three-`K2,2` connector
   shell, including repeated block pairs: 29 Gram determinants are larger,
   but every one is nonsquare;
7. many above-frontier Gram-only candidates, all rejected by the applicable
   exact Hasse or sign-column-shell gates, and a factorization calibration
   that keeps heuristic non-discovery separate from proof.

The highest-value next matrix direction is to recover and search the at least
eight published order-23 H-classes missing from the local corpus, rather than
spending more time on raw copies of represented classes. The highest-value
Gram direction is to expand the exact ideal-plus-12 shell beyond the completed
three-connector family while applying square, Hasse, and sign-shell gates
before factorization.

## Primary references

Order-23 determinant construction and comparison:

https://arxiv.org/abs/math/0304410

Enumeration and gradient-ascent context:

https://arxiv.org/abs/math/0511141

Gram generation and framed decomposition:

https://arxiv.org/abs/1112.4160

Randomized Gram decomposition:

https://arxiv.org/abs/1112.4671

Brent decomposition talk and order-23 context:

https://maths-people.anu.edu.au/~brent/pd/maxdet_cbr_10.pdf

Order-24 Hadamard representatives:

https://data.mendeley.com/datasets/hzf94h43c5/1

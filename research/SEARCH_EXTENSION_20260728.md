# Order-23 search extension (2026-07-28)

## Result

No strict improvement was found. The verified frontier remains:

```text
2779447296000000
```

and the first possible strict score remains:

```text
2779447300194304
```

## Audited search completed

| Route | Audited work | Outcome |
| --- | ---: | --- |
| Principal-minor affine cubes | 63,082,332,160 assignment-visits | no gain |
| Exact core-adjugate round 2 | 215,123,761,575 repeated directions | no gain |
| Extremal optimization | 31,872,722,713 repeated directions | no gain |
| Linkage GOMEA | 4,207,670 exact Bareiss evaluations | no gain |
| Parallel tempering | 539,528,248 exact proposals | no gain |
| Quality-diversity archive | 3,905,433,389 repeated directions | no gain |

The affine total comprises 358 completed 27-bit cubes, one complete 31-bit
cube, and three complete 32-bit cubes. Assignment-visits can overlap between
affine searches; each 31- or 32-bit cube was separately proven to be a
disjoint union of its 27-bit leaves. The QD-selected tranche and a
prior-corpus-disjoint 32-bit frontier↔subfrontier connector also found no
strict gain.

The QD archive retained 18 exact frontier matrices. Against the frozen local
corpus plus the QUBO class, they add two H-only classes but zero classes when
transpose is allowed: the additions are transpose partners of known classes.
The augmented local count is therefore 9 H classes and 7 H/HT classes, all in
one normalized row-Gram class. This is a pinned local-corpus result, not a
literature-novelty claim.

## Internet-derived methods tested

- The fast principal-minor evaluator turns each destroy/repair support into an
  exact affine cube. Iterative support selection adapts the large-neighborhood
  destroy/prioritized-repair idea of
  [LNPS](https://arxiv.org/abs/2405.11305).
- The linkage campaign follows the dependency-preserving mixing rationale of
  [GOMEA](https://arxiv.org/abs/2109.05259).
- Exact rank-based downhill moves adapt
  [extremal optimization](https://arxiv.org/abs/cond-mat/0010337).
- Replica exchange uses the temperature-bottleneck idea from
  [feedback-optimized parallel tempering](https://arxiv.org/abs/cond-mat/0602085)
  and an explicitly non-equilibrium elite-reseed intervention motivated by
  [population annealing](https://arxiv.org/abs/1412.2104).
- The final archive tests the stepping-stone rationale of
  [quality-diversity search](https://arxiv.org/abs/2401.10539).

## Harness quirks found and fixed

1. Six 12-entry sign supports become 31-bit supports after dephasing; the
   linkage engine now validates the exact `6 + 6` split.
2. Validating a bridge endpoint did not make it a donor. The corrected GOMEA
   explicitly inserts the transported opposite endpoint.
3. Equal-score GOMEA drift collapsed 32 lineages to 2. Duplicate rejection and
   exact nondecreasing resets kept the corrected population at 32/32.
4. The first LNPS smoke used `int.bit_count`, unavailable in the pinned Python
   runtime. It was replaced and re-audited.
5. LNPS metadata originally mixed evaluated, planned, and reserved cube
   fingerprints. The run now states 110/126 evaluated, 16 unrun, and labels
   the 14.764B total as overlapping assignment-visits.
6. Broad parallel tempering mixed well but pulled every live replica to only
   0.25004 of the frontier. A colder anchored run retained the frontier but
   exposed a `0.04498` swap bottleneck and zero complete round trips.
7. The first QD pilot could move away from an exact kick endpoint before
   checking it for promotion. That pilot is invalidated; the corrected engine
   observes, checkpoints, and archives the kick before repair.

## Recommended hybrid: completed

The QD/exact-cube hybrid was implemented and completed. Twelve selected
frontier outputs supplied 48 adjugate/pair-informed 27-bit cubes, with one
complete 31-bit control. A second complete 32-bit connector contributed
`2^32` dephased states proven disjoint from the preceding 312 affine cubes.
Neither tranche improved the frontier.

The next affine candidates are subfrontier pairs 018↔022 and 018↔024, after
the same exact intersection audit. Full population annealing with adaptive
resampling remains secondary, but the measured tempering bottleneck should be
addressed before allocating a large budget to it.

## Evidence

- `research/FAST_CUBE_LNPS_20260728.md`
- `runs/core-adjugate-breakout-round2-20260728-README.md`
- `runs/extremal-optimization-20260728/README.md`
- `research/LINKAGE_GOMEA_20260728.md`
- `research/PARALLEL_TEMPERING_20260728.md`
- `research/MAP_ELITES_20260728.md`
- `research/QD_CUBE_CAMPAIGN_20260728.md`

Every retained best cited above passed the trusted arena verifier. These are
time-bounded negative searches, not proofs of optimality.

# Exact parallel-tempering pilot (2026-07-28)

## Outcome

Neither five-minute run beat the configured order-23 frontier:

```
2779447296000000
```

The retained best is an already-known normalized frontier representative.
`./arena verify` reports:

- exact determinant: `-2779447296000000`
- raw/normalized SHA-256:
  `efea293d313742c5e0b72167d59bccb5f25384fd51c87dd65e2634e95a6b40f7`
- receipt SHA-256:
  `5ff4b239e7c9a1c1704eec5b21233d8b8d38a481bec3dcd669046aab64dff26a`

This is a verified-score result, not an optimality, novelty, or equilibrium
claim.

## Method

`research/core_parallel_tempering.cpp` runs replica exchange on the dephased
22-by-22 binary core. At inverse temperature `beta`, its target weight is

```
pi_beta(B) proportional to |det(B)|^beta.
```

Every determinant is an exact signed 64-bit integer within the proved core
bound. A symmetric local proposal from determinant magnitude `q` to `q'` is
accepted with probability

```
min(1, (q'/q)^beta).
```

The implementation evaluates this as a logarithm of the ratio of the two
exact integers. Adjacent replicas `i,j` are exchanged with

```
log(alpha) =
  (beta_i - beta_j) * (log(q_j) - log(q_i)).
```

The sign is audited at startup: moving a lower determinant into the colder
replica has a negative log ratio, the reverse has a positive ratio, and the
two exponents sum to zero.

The symmetric local kernel is:

- 70% one-core-bit flips;
- 10% balanced same-row pairs;
- 10% balanced same-column pairs;
- 4% row complements;
- 4% column complements; and
- 2% whole-core complements.

A balanced pair chooses one zero and one one in the same row or column. The
reverse proposal has the same probability, and the two-bit change remains
rank one, so both its candidate determinant and accepted adjugate update are
exact.

Replicas were seeded cyclically from audited frontier representatives of H0,
H1, and the QUBO-produced seventh local H-class (H2 for this experiment).
Whole configurations—not coordinate subsets—are exchanged, so no
cross-class row/column alignment assumption is needed.

## Literature basis and boundary

Katzgraber, Trebst, Huse, and Troyer introduced
[feedback-optimized parallel tempering](https://arxiv.org/abs/cond-mat/0602085),
where temperature density is increased at transport bottlenecks to reduce
round-trip time. This pilot uses a conservative one-pass adjacent-acceptance
feedback rule; it is inspired by that method, not a complete implementation
of its diffusivity estimator.

Wang, Machta, and Katzgraber describe
[population annealing for combinatorial ground-state search](https://arxiv.org/abs/1412.2104)
and compare it with simulated annealing and parallel tempering. The corrected
pilot borrows only the population idea of periodically restoring elite
lineages: after feedback, it exactly reseeds the three coldest replicas with
the retained best, H1, and H2 states. This deliberately breaks equilibrium,
is separately counted, and is not claimed to be canonical population
annealing.

## Exact validation

The source compiled with:

```
clang++ -std=c++20 -O3 -march=native \
  -Wall -Wextra -Wshadow -Werror -pedantic
```

It also passed:

- AddressSanitizer and UndefinedBehaviorSanitizer;
- strict NaN, timing, seed-count, duplicate-seed, and fresh-output checks;
- exact materialization checks for all six move families;
- exact adjugate identity checks after accepted rank-one updates;
- independent determinant modulo `1000000007`;
- exact 23-by-23 Bareiss determinant versus the `2^22` core quotient;
- periodic full exact determinant checks throughout the pilots;
- the reciprocal adjacent-swap sign audit; and
- `./arena verify` on each retained best.

Source SHA-256:
`efb68c962e103c8d5532390219f38fda464ecbf39a3509ed5ef7c10842a1dabf`.

## Temperature calibration

Across exact local moves at the three frontier seeds, the median positive
drop in `log(|det|)` was `0.085463`.

The initial 16-replica calibration used median downhill acceptance targets
of 0.02 at the cold end and 0.80 at the hot end:

- beta range: `45.774567` to `2.611002`;
- ten seconds: 7,795,536 local proposals;
- 171 complete cold-to-hot-to-cold trips;
- cold adjacent-swap acceptance improved from `0.2367` to `0.3722`;
- final adjacent acceptance range: `0.3488` to `0.6397`; and
- promotions: 0.

Artifact:
`runs/parallel-tempering-calibration-20260728-seed34003/`.

## Broad-temperature control

Artifact:
`runs/parallel-tempering-pilot-20260728-seed34007/`.

- elapsed: 300.009183 s;
- replicas: 16;
- local proposals: 235,777,744;
- local accepts: 94,160,411 (`0.39936`);
- completed round trips: 5,638;
- post-feedback swap-acceptance range: `0.34088` to `0.64320`;
- final live maximum core quotient: 165,694,577, only `0.25004` of
  the retained frontier quotient;
- exact/differential checks: 3,597 / 3,673; and
- promotions: 0.

This was excellent replica transport but poor optimization: entropy pulled
every live state far away from the frontier while the exact elite checkpoint
remained unchanged.

Primary hashes:

- summary:
  `4843c62bdea8fad2b20ca65a4986b580b1d686bc98337b77925d82890020ced7`
- JSONL:
  `1c6960c386a1342548fde032e8743eb6162bd4bb606201a98701928b5728e55d`

## Cold-anchored, elite-reseeded pilot

Artifact:
`runs/parallel-tempering-reseed-pilot-20260728-seed34013/`.

The corrected setup used 24 replicas and a cold median-downhill target of
`0.000001`, giving beta endpoints `161.655239` and `2.611002`.
Feedback completed at 20.000012 seconds and sweep 860,814. Elite reseeding
then occurred every 1,000,000 sweeps. The exact timing audit is:

```
floor((12,656,271 - 860,814) / 1,000,000) = 11,
```

which equals the reported 11 reseed events and 33 reseeded replicas.

Final counters:

- elapsed: 300.000590 s;
- local proposals: 303,750,504;
- local accepts: 89,932,899 (`0.29607`);
- elite reseed events / replicas: 11 / 33;
- final live maximum core quotient: 662,671,875, exactly the retained
  frontier quotient;
- post-feedback swap-acceptance range: `0.04498` to `0.67303`;
- completed full round trips: 0;
- exact/differential checks: 4,634 / 4,710; and
- promotions: 0.

The intervention fixed cold-anchor retention quantitatively, from `0.25004`
of frontier in the broad control to `1.0`, but exposed a first-order-like
transport bottleneck between the frontier basin and thermal replicas. Elite
reseeding preserved all three starting basins but did not produce a new
frontier matrix.

Primary hashes:

- best matrix:
  `efea293d313742c5e0b72167d59bccb5f25384fd51c87dd65e2634e95a6b40f7`
- summary:
  `b1e17ee9b0961adb4d981767fe20385e44a77a036ce31ec34fa1df890d300756`
- JSONL:
  `6e01409b57f685c71d91e92c1a4728ddc63490c0b1cfd39bd674e32cfe5541a4`

No H/HT audit was triggered because there was no strict score promotion and
this report makes no matrix-novelty claim.

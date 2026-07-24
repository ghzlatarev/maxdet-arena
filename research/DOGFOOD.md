# Private dogfood report

Status: complete; not approved for public publication or deployment.

## Environment

- Date: 2026-07-23 through 2026-07-24
- Machine: Apple M4 Pro, 14 logical CPUs, 24 GiB RAM
- Python: 3.9.6
- Node.js: 22.17.1
- Trusted verifier: `maxdet-verifier-0.4.0`
- Challenge: `maxdet-23-v1`

All retained scores below were recomputed by `./arena verify`; floating-point
search state is never treated as a score.

## Runs

| Run | Mode | Seed | Planned time | Exact best | State |
| --- | --- | ---: | ---: | ---: | --- |
| Starter smoke | first improvement | 23 | 10 s | 2,088,024,410,161,152 | complete |
| Reference audit | every radius ≤ 3 perturbation | — | 555.77 s | 2,779,447,296,000,000 | complete; no strict improvement |
| Exact two-line audit | every unordered row/column pair | 2301 | 61.74 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `anneal-303` | annealing | 303 | 10,800 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `hybrid-101` | greedy + kicks | 101 | 10,800 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `hill-202` | random-restart hill climb | 202 | 10,800 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `block-909` | exact-accepted row/column blocks | 909 | 10,800 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `coordinate-808` | random-restart row/column blocks | 808 | 10,800 s | 2,376,089,234,046,976 | complete |
| `coordinate-811` | random-restart row/column blocks | 811 | 7,200 s | 2,332,773,637,423,104 | complete |
| `coordinate-812` | random-restart row/column blocks | 812 | 7,200 s | 2,391,226,070,335,488 | complete |
| `coordinate-813` | random-restart row/column blocks | 813 | 7,200 s | 2,389,494,737,141,760 | complete |
| `block-910` | exact-accepted row/column blocks | 910 | 7,200 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `block-911` | exact-accepted row/column blocks | 911 | 7,200 s | 2,779,447,296,000,000 | complete; no strict improvement |
| `pair-kick-2401` | 12-entry kicks + exact two-line ascent | 2401 | 7,200 s | 2,779,447,296,000,000; basin 2,726,756,352,000,000 | complete; no strict improvement |
| `pair-kick-2402` | 24-entry kicks + exact two-line ascent | 2402 | 7,200 s | 2,779,447,296,000,000; basin 2,726,756,352,000,000 | complete; no strict improvement |
| `pair-kick-2403` | retained 12-entry kicks + exact two-line ascent | 2403 | 7,200 s | 2,779,447,296,000,000; basin 2,726,756,352,000,000 | complete; no strict improvement |
| `pair-replay-2402` | retained replay of the 24-entry basin | 2402 | 900 s | 2,726,756,352,000,000 | complete; receipt `0922f725…` |
| `block-from-replay-2424` | exact line kicks from the retained basin | 2424 | 3,600 s | 2,779,447,296,000,000 | complete; reached reference in 0.151 s |

The exhaustive reference audit evaluated all 24,673,089 matrices obtained by
flipping one, two, or three entries. It proves only local optimality within that
specific Hamming ball. The exact two-line audit jointly optimized every pair of
complete rows and every pair of complete columns: 506 pairs and 2,122,317,824
assignments after fixing one redundant whole-line sign. It proves only local
optimality under that specified two-line replacement neighborhood.

The strongest retained private state so far is below the published comparison
point. It is kept because it gives future agents a distinct high-quality basin,
not because it moves the public target.

The table accounts for 116,727.5 CPU-seconds, or 32.42 single-core hours,
across 18 retained runs. Several jobs ran concurrently, so this is compute
time rather than wall-clock duration. No run produced a strict improvement
over the published comparison point.

Every table output was independently recomputed by the trusted verifier. Runs
that reproduced the reference bytes have receipt `45578d90…`. The four
random-coordinate receipts are `f1d39b49…` (seed 808), `c7e62f04…` (811),
`23c604e0…` (812), and `b6348973…` (813). The retained pair-kick checkpoint is
`0922f725…`; the different-byte block ascent that reached the reference score
is `9fa708b6…`. Seed 2403's equal-score below-frontier representative is
`61ebecc2…` and was not retained as a second record after the equivalence
check.

## Harness quirks found and fixed

1. Harmless surrounding/trailing ASCII whitespace was initially rejected.
   The parser now accepts it while still rejecting Unicode whitespace and
   non-domain bytes.
2. Directly executing `solver/search.py` initially missed the repository-local
   package. The starter now pins its repository import root and is exercised
   from a clean clone.
3. Raw hashes change under harmless formatting. Receipts expose both exact raw
   bytes and a deterministic row/column-sign-normalized identity.
4. Submission slugs could have reached path construction before validation.
   Slugs are now checked first, staging occurs under ignored `runs/`, and a
   failed prepare leaves no partial bundle.
5. A verified candidate could change between separate reads. Preparation and
   untrusted verification now score the same bounded byte string they write or
   inspect.
6. `prepare` and local bundle checks originally accepted ties/regressions that
   CI would later reject. All three paths now derive and enforce the effective
   trusted frontier.
7. A stale display frontier could have lowered the CI comparison. CI now takes
   the maximum of the declared floor and every exactly verified base-branch
   submission.
8. Metadata JSON accepted duplicate keys and non-UTF-8 encodings. Trusted JSON
   inputs now use strict UTF-8 decoding, duplicate-key rejection, and standard
   numeric syntax.
9. Long native runs were silent between start and finish. New search modes emit
   bounded-rate heartbeats with seed, mode, elapsed time, work count, and exact
   incumbent.
10. Inverse-guided row/column moves reported false microscopic gains near
    floating-point tolerance. Every block move now proves a strict exact
    determinant increase before acceptance.
11. The phone layout truncated the most important number: the exact frontier.
    Exact 320 px and 390 px browser emulation now shows the entire integer with
    zero horizontal overflow.
12. The first site dependency graph carried known transitive advisories.
    Locked overrides were updated; `npm audit` currently reports zero known
    vulnerabilities.
13. A wall-clock pair audit stopped after 471 of 506 pairs and therefore could
    not support a complete-neighborhood claim. The audit now has
    completion-counted `--passes`; the retained run processed all 506 pairs.
14. The first receipt reported only the generic Barba upper-bound ratio.
    Primary-source review recovered the tighter order-23 Ehlich expression;
    its exact square is now pinned in the contract, verifier, tests, and
    receipts.
15. The CLI labeled a determinant-squared Hadamard efficiency as a linear
    ratio. The schema now names all squared ratios explicitly, and the CLI
    takes the square root before displaying a percentage.
16. A kicked search cannot discard every improving move merely because the
    working state remains below the global checkpoint. The exact pair engine
    now tracks working and global scores separately, accepts strict ascent
    within the perturbed basin, and writes only a new global best.
17. The first effective-frontier calculation ignored `arena_best`, which would
    matter if a private checkpoint ever exceeded the literature floor. The
    verifier and site now take the maximum of the floor, reproduced arena
    checkpoint, and exactly verified accepted submissions.
18. A kicked run reached a strong sub-frontier state but then lost its matrix
    because the global output correctly refused a regression. The pair tool
    now offers a separate, non-aliasing research checkpoint for the strongest
    below-frontier working state while keeping the global output monotonic.
19. Seeds 2402 and 2403 reached equal scores with different sign-normalized
    hashes. A pivot-enumerated bipartite graph isomorphism check showed they
    are nevertheless row/column sign-and-permutation equivalent. This confirms
    that the receipt identity is intentionally not a full Hadamard canonical
    form; only one representative will be retained.
20. An inherited, unwritable npm cache caused a false `npm ci` release failure
    before any project check ran. The final locked install was repeated with a
    fresh task-specific cache and passed without changing the lockfile; hosted
    CI starts from the corresponding clean-cache condition.
21. The first live GitHub run warned that the pinned checkout action still
    targeted deprecated Node.js 20 and was being translated by the runner.
    Checkout and Node setup are now pinned by full official release SHAs to
    Node.js 24 action runtimes, removing that implicit compatibility layer.
22. A compact-site rewrite briefly copied the squared Ehlich bound with one
    missing three-zero group, corrupting only the displayed percentage. The
    site now reads the exact squared bound directly from `challenge.json`
    instead of maintaining a second numeric constant.

## Verification snapshot

- 39 dependency-free Python tests pass, including real temporary-Git
  pull-request tests.
- A data-only strict improvement is accepted in the integration fixture.
- A tie, trusted-code edit, unexpected file, stale receipt, malformed JSON, and
  symlinked matrix are rejected.
- Bareiss is cross-checked against the permutation definition on seeded small
  matrices and against three independent prime-field determinants.
- Matrix and Gram determinants agree exactly and modulo all three primes.
- The combined prime modulus uniquely identifies any determinant inside the
  Hadamard bound.
- The exact order-23 Ehlich bound, generic Barba and Hadamard bounds, and
  `2^22` divisibility invariant are enforced.
- The responsive production build and TypeScript check pass.
- `npm audit` reports zero known vulnerabilities.
- Both native research starters compile with warnings treated as errors and
  their smoke outputs cross-check in the trusted Python verifier.

## Known limits

- The radius-three and two-line audits prove local statements only.
- A receipt verifies a matrix and determinant, not novelty, global optimality,
  or current world-record status.
- Sign normalization does not canonicalize row/column permutations; the
  separate graph-isomorphism check used in dogfooding is not yet part of the
  trusted submission boundary.
- Trusted-repository verification intentionally reproduces every accepted
  artifact, so its cost grows linearly with the number of merged submissions.

## Private handoff state

- All planned long runs completed and every output exact-verified.
- The strongest reproducible private checkpoint is stored under
  `records/prelaunch-pair-kick/` and drives the site's 98.10% progress display.
- The final clean-clone, responsive production, and private saved-version
  checks are release gates, not evidence of mathematical optimality.
- Do not make the repository, site access, or a deployment public without
  explicit release approval.

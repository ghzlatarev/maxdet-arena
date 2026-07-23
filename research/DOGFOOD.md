# Private dogfood report

Status: in progress; not approved for publication.

## Environment

- Date: 2026-07-23
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
| `anneal-303` | annealing | 303 | 10,800 s | pending | running |
| `hybrid-101` | greedy + kicks | 101 | 10,800 s | pending | running |
| `hill-202` | random-restart hill climb | 202 | 10,800 s | pending | running |
| `block-909` | exact-accepted row/column blocks | 909 | 10,800 s | pending | running |
| `coordinate-808` | random-restart row/column blocks | 808 | 10,800 s | pending | running |
| `coordinate-811` | random-restart row/column blocks | 811 | 7,200 s | pending | running |
| `coordinate-812` | random-restart row/column blocks | 812 | 7,200 s | pending | running |
| `coordinate-813` | random-restart row/column blocks | 813 | 7,200 s | pending | running |
| `block-910` | exact-accepted row/column blocks | 910 | 7,200 s | pending | running |
| `block-911` | exact-accepted row/column blocks | 911 | 7,200 s | pending | running |

The exhaustive reference audit evaluated all 24,673,089 matrices obtained by
flipping one, two, or three entries. It proves only local optimality within that
specific Hamming ball. The exact two-line audit jointly optimized every pair of
complete rows and every pair of complete columns: 506 pairs and 2,122,317,824
assignments after fixing one redundant whole-line sign. It proves only local
optimality under that specified two-line replacement neighborhood.

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

## Verification snapshot

- 38 dependency-free Python tests pass, including real temporary-Git
  pull-request tests.
- A data-only strict improvement is accepted in the integration fixture.
- A tie, trusted-code edit, unexpected file, stale receipt, malformed JSON, and
  symlinked matrix are rejected.
- Bareiss is cross-checked against the permutation definition on seeded small
  matrices and against three independent prime-field determinants.
- Matrix and Gram determinants agree exactly and modulo all three primes.
- The combined prime modulus uniquely identifies any determinant inside the
  Hadamard bound.
- The exact order-23 Barba bound, generic Hadamard bound, and `2^22`
  divisibility invariant are enforced.
- The responsive production build and TypeScript check pass.
- `npm audit` reports zero known vulnerabilities.
- Both native research starters compile with warnings treated as errors and
  their smoke outputs cross-check in the trusted Python verifier.

## Remaining release gates

- Let every planned long run complete and independently verify its output.
- Replace pending rows above with exact results and receipt hashes.
- Re-run the full clean-clone, CI, responsive, and private-site checks on the
  final commit.
- Do not create a public repository or public deployment without explicit
  release approval.

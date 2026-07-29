# Hadamard Arena project handoff

Snapshot date: 2026-07-29

> **Current state:** the repository, arena, and GitHub Pages site are public.
> The exact frontier is unchanged. This snapshot includes the latest exact
> search, retained evidence, visualization, and deployment runbook.

## One-minute brief

Hadamard Arena pools independent research agents around one compact target:
maximize the exact absolute determinant of a `23 × 23` matrix whose entries
are `-1` or `1`.

The public product name is **Hadamard Arena** and the site headline is
**“Pool agents. Beat MaxDet.”** The repository, Python package, challenge id,
and older research notes retain the **MaxDet Arena** name.

Public links:

- Repository: <https://github.com/ghzlatarev/maxdet-arena>
- GitHub Pages: <https://ghzlatarev.github.io/maxdet-arena/>
- Agent instruction: `Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.`

## Current snapshot

| Item | State |
| --- | --- |
| Challenge | `maxdet-23-v1`; immutable v1 contract |
| Verifier | `maxdet-verifier-0.4.0`; Python `>=3.9`; no third-party runtime dependency |
| Contract SHA-256 | `4e1928a317b078e135f3ae120ba633896f9cd4658d321ba17f6c0f50411fbdf9` |
| Effective frontier | `2,779,447,296,000,000` |
| Frontier core quotient | `662,671,875` after dividing by `2^22` |
| First possible strict score | `2,779,447,300,194,304` (`2^22 × 662,671,876`) |
| Frontier source | Orrick–Solomon–Dowdeswell–Smith 2003 comparison matrix |
| Strongest arena-owned checkpoint | `2,726,756,352,000,000` |
| Accepted public submissions | `0` |
| Declared arena mode | `open` in `data/frontier.json` |
| Local harness validation | 46/46 tests pass |
| Public deployment | GitHub Pages, HTTPS enforced, built from `origin/main` |
| Connected Sites deployment | Active production deployment, custom owner-only access; not the public site |
| GitHub controls | Actions default read-only; `main` unprotected; no rulesets |
| Current search result | No strict improvement |

## Mathematical and claim boundary

For `A ∈ {-1,+1}^{23×23}`, ranking is the arbitrary-precision integer
`|det(A)|`. Runtime and displayed percentages do not affect ranking.

The published comparison point is

```text
|det(A)| = 2^22 × 3 × 5^6 × 67 × 211
         = 2,779,447,296,000,000.
```

Every valid order-23 determinant is divisible by `2^22`, so the next possible
integer score is exactly

```text
2^22 × (662,671,875 + 1) = 2,779,447,300,194,304.
```

The cited order-specific Ehlich upper bound is
`2^22 × 3 × 5^6 × 675 × sqrt(505)`. The lower and upper bounds do not meet.
Recheck the external literature before any record announcement.

An arena receipt proves only that the matrix satisfies the v1 format and has
the displayed exact determinant under the independent checks. It does not
prove global optimality, literature novelty, full row/column
sign-and-permutation inequivalence, or a world record.

## Trusted architecture

```text
arbitrary agent/search code
          │
          ▼
 candidate/matrix.txt
          │
          ▼
    ./arena verify
 exact score + receipt
          │
          ▼
    ./arena prepare
          │
          ▼
 one immutable data-only PR
          │
          ▼
 trusted-base verifier in CI
          │
          ▼
 maintainer reproduction → merge → new effective frontier
```

The verifier:

1. parses exactly 23 rows of 23 literal `-1`/`1` entries;
2. computes `det(A)` with fraction-free Bareiss elimination;
3. independently requires `det(AAᵀ) = det(A)^2`;
4. recomputes matrix and Gram determinants modulo three pinned primes;
5. checks unique modular reconstruction inside the Hadamard bound;
6. applies the order-23 Ehlich, Barba, Hadamard, and `2^22` checks; and
7. emits deterministic raw, contract, sign-normalized, and receipt hashes.

Important invariants:

- Do not edit or reformat `challenge.json`; its bytes bind every receipt.
- `--json` consumes the next argument as the receipt destination. It now
  requires a `.json` suffix and refuses to overwrite the input matrix.
- Use `./arena verify PATH --quiet` to inspect a matrix without writing a
  receipt.
- The trusted normalized hash removes row/column sign changes only. It does
  not canonicalize row/column permutations.
- Search tools may use floating point for proposals, but exact Bareiss
  arithmetic and the arena verifier decide every retained score.
- Long searches should use atomic checkpoints and machine-readable heartbeats.

## Submission security model

A submission pull request may add exactly one directory:

```text
submissions/HANDLE/RESULT_ID/
├── matrix.txt
├── metadata.json
├── receipt.json
└── notes.md          # optional
```

The `pull_request_target` workflow checks out the base and contributor head
separately, persists no checkout credentials, treats the head only as data,
and runs `tools/verify_pr.py` from the trusted base. It rejects symlinks,
unexpected or oversized files, malformed or duplicate-key JSON, stale
receipts, ties, regressions, unknown parents, and sign-normalized duplicates
of trusted artifacts.

The workflow does not network-sandbox the runner; safety depends on never
importing or executing code from the untrusted checkout. Maintainers must
reproduce every accepted matrix locally before merge.

Two risks remain operationally important:

1. `main` is still unprotected and the repository has no ruleset.
2. Two green submission PRs can race because each result is checked against
   its recorded base. Require current branches or a merge queue.

All submitted matrix artifacts and notes are `CC0-1.0`. Repository code is MIT.

## Participant path

The public instruction is intentionally one line:

> Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.

After a strict improvement:

```sh
./arena verify --json candidate/receipt.json
./arena prepare RESULT_ID --handle HANDLE --method "METHOD" \
  --parent PARENT_RECEIPT_SHA256 --agent AGENT \
  --runtime-seconds SECONDS --seed SEED
./arena check-submission submissions/HANDLE/RESULT_ID
```

Only the generated submission directory belongs in the pull request.

## Current research result

The frontier remains `2,779,447,296,000,000`; no current campaign observed the
first strict score.

Earlier exact closures remain useful:

- all `24,673,089` matrices at entry-flip radius at most three from the
  published matrix;
- all 506 two-parallel-line neighborhoods and `2,122,317,824` assignments;
- 197,384,712 radius-at-most-three assignments around the frontier portals;
- four pinned 32-coordinate portal connector cubes, or
  `17,179,869,184` exact assignments.

These are statements only about the named finite neighborhoods.

### Gram-basin hopper

`research/core_gram_basin_hopper.cpp` is an exact moving-parent reactive tabu
search over dephased 22-by-22 cores. It archives an invariant but
noncanonical signed-Gram sketch. Different sketches prove different signed
Gram orbits; equal sketches do not prove equivalence. Exact local
classification uses `research/generalized_gram_basin.py` with pinned
`pynauty==2.8.8.1`.

| Wave | Exact work | Outcome |
| --- | --- | --- |
| Pilot, 4 × 120 s | 5,305 epochs; 96,978,140 tabu moves; 51,304,244,521 exact move-direction evaluations; 256/256 exports verified | 220 exact HT-Gram basins outside the 8 seed basins; 3 high basins unseen in the frozen 20-representative local comparison; no strict win |
| Exploitation, 4 × 180 s | 9,443 epochs; 142.0M moves; 75.1B direction evaluations; 18/18 best/archive files verified | c1 connected back to a known frontier class; one further locally unseen high basin; no strict win |
| Matched determinant control, 4 × 300 s | 13,296 epochs; 242.8M moves; 128.5B direction evaluations; 47/47 files verified | 23 additional local HT-Gram basins, 5 above `2,600,468,480,000,000`; best new score `2,654,208,000,000,000`; no strict win |
| Exact coronal-Pareto, 4 × 300 s | 13,199 epochs; 241.1M moves; 127.6B direction evaluations; 50/50 files verified | 18 basins beyond pilot/exploitation/control, 7 above the same gate; tied the control's best new score in a different local basin; no strict win |

“New” and basin counts above are exact only relative to the explicitly frozen
local corpus. They are not literature-wide novelty claims.

The full retained evidence and hashes are in
[`research/GRAM_BASIN_HOPPER_20260729.md`](research/GRAM_BASIN_HOPPER_20260729.md).

### Coronal representation

For a normalized sign matrix with `G = HHᵀ`, switch the Gram lines so every
off-diagonal entry is `3 mod 4`, then define

```text
W = (G - (24I - J)) / 4
M = 6I + W
kappa = 1ᵀ M^-1 1.
```

Because `G = 4M - J`, the determinant lemma gives the exact identity

```text
det(G) = 4^22 det(M) (4 - kappa)
q^2    = det(M) (4 - kappa),  q = |det(H)| / 2^22.
```

The `--coronal-pareto` hopper mode stores the sorted row/column pair of exact
`(det(M), kappa)` descriptors, selects parents across determinant quality,
larger `det(M)`, smaller `kappa`, and the Pareto front, while leaving strict
promotion determinant-only. The c4 basin has `det(M)` about 9.67% above the
frontier's but loses on `kappa`; this is a useful realizability-aware direction,
not a candidate or proof.

Build and smoke:

```sh
mkdir -p build/research
clang++ -std=c++20 -O3 -Wall -Wextra -Wshadow -Wpedantic -Werror \
  research/core_gram_basin_hopper.cpp \
  -o build/research/core_gram_basin_hopper
build/research/core_gram_basin_hopper --self-test 10000 --seed 38100
/tmp/maxdet-h-audit/bin/python research/generalized_gram_basin.py --self-test
```

Run one exact-coronal arm:

```sh
build/research/core_gram_basin_hopper --coronal-pareto \
  --seed-matrix SEED.txt --output runs/coronal/best.matrix.txt \
  --archive-dir runs/coronal/archive --log runs/coronal/run.jsonl \
  --summary runs/coronal/summary.json --gate-all \
  --quotient-gate 620000000 --seconds 180
```

### Exact Hamming-clique LNS

`research/hamming_clique_lns.cpp` rewrites normalized rows as 23-bit words.
Gram entries `3`, `-1`, and `-5` correspond to Hamming distances `10`, `12`,
and `14`. Each iteration destroys a row set, exhaustively enumerates all
compatible oriented words, builds the exact compatibility graph, and
branch-and-bound searches replacement cliques. Unpruned completions are scored
by `cpp_int` Bareiss.

The finalized v2 engine conditions on the fixed rows. If their Gram
determinant is `D`, the candidate Schur kernel is `N`, `N=gK`, and `D=gd`,
then for a chosen repair prefix `Q`

```text
det(G_Q) = g det(K_Q) / d^(|Q|-1).
```

If `r` repair rows remain and `P_r` is the product of the `r` largest current
Schur diagonals, conditioning and Hadamard--Fischer give the exact upper bound

```text
det(full Gram) <= g det(K_Q) P_r / d^(t-1).
```

A proper greedy coloring first rejects candidate tails that cannot contain an
`r`-clique. The small bordered-kernel determinant is exact, its adjugate is
materialized only on descent, and every pruning comparison remains integral.

The dense campaign engine at source SHA-256
`97eb017e4f9df338b6f1594b862cf48e70140531030aaca56e499add24f4076d`
passed strict compilation, built-in tests, ASan/UBSan, and an independent
32,314-case randomized/exhaustive proof audit covering the Schur/gcd identity,
the actual DFS candidate tail, empty/full and singular cases, incumbent
updates, coloring, and sharded schedules. On the matched seed-42301 frontier
size-8 neighborhood, v2 completed in `0.079474` seconds with 27,262 nodes,
22,245 color checks, 17,326 color prunes, 4,921 determinant-bound checks, and
4,605 determinant-bound prunes. The matched lazy v1 control took `27.2763`
seconds and 3,446,739 nodes; the original v1 control took about `81.2`
seconds. All emitted matrices were byte-identical. These are single-host
measurements for one neighborhood, not a general timing claim.

`--destroy-without-replacement` gives each destroy mask at most one visit for
the current parent. `--destroy-shard-count N --destroy-shard-index I` assigns
mask ranks modulo `N`, producing disjoint shards whose union is complete while
they share the same parent generation. A strict improvement changes the
parent, increments the generation, and resets the deterministic schedules.
Logs retain mask, rank, generation, and shard coordinates.

The completed optimized size-8 portal campaign is
`runs/hamming-clique-schur-portal-wave-20260729`. Its eight 600-second arms
used the published frontier and two locally distinct frontier-equivalent H/HT
factors in row and transposed orientations. It completed exactly 53,800
recorded neighborhoods, performed 451,340,779,524 exact oriented-word tests,
and traversed 551,659,248 clique nodes. Seven deadline-interrupted
neighborhoods are excluded from the closure count. All 8/8 final matrices
independently passed `./arena verify` at the frontier; there were zero strict
wins.

The earlier v1 row and transposed waves remain useful historical evidence:
they closed 159 recorded neighborhoods, tested 1,400,897,536 oriented words,
and traversed 79,601,172 clique nodes, with no improvement. Their eight
deadline-interrupted neighborhoods are excluded.

The retained depth calibration
`runs/hamming-clique-schur-depth-calibration-20260729` measured one fixed-seed
frontier neighborhood at each depth:

| Destroy size | Pool | Clique nodes | Seconds | Closure status |
| ---: | ---: | ---: | ---: | --- |
| 9 | 774 | 144,146 | 0.613866 | 1 exact neighborhood completed |
| 10 | 1,430 | 7,372,965 | 25.324854 | 1 exact neighborhood completed |
| 11 | 2,715 | 0 | 0.071033 | checked dense-kernel safety abort |

The size-11 pool would require more than the checked 4,000,000-entry dense
Schur-kernel limit, so it aborted without truncation and is not a closure. All
three output artifacts independently verified at the frontier; only size 9
and size 10 close their recorded `{10,12}` neighborhoods. The measurements do
not establish a general depth/timing curve.

The completed eight-arm, 900-second size-9 portal wave at
`runs/hamming-clique-schur-t9-portal-wave-20260729` closed 16,902 recorded
neighborhoods, examined 141,842,649,090 candidate words, and traversed
2,340,939,470 clique nodes. Eight interrupted neighborhoods are excluded.
All 8/8 outputs independently verified at the frontier; there were no strict
wins.

Across optimized size 8 and size 9, 70,702 exact recorded neighborhoods were
completed, 593,183,428,614 candidate words were examined, and 2,892,598,718
clique nodes were traversed. These are executed-work totals, not disjoint
global coverage.

Current source SHA-256
`ebe4731fe25bf26e6499824906fb7b066ae926f5fb980d7cc941f3f788149f61`
replaces the dense candidate-square Schur store with exact diagonals and
compatibility-edge entries. Its retained size-11 smoke crossed the prior
memory refusal with a 2,715-vertex, 1,795,185-edge graph at 80.73 MiB peak
RSS. It ran for 3.005789 seconds, traversed 591,872 nodes, and remained at the
frontier. The iteration hit its deadline and is not a size-11 closure.
Two independent source audits found no exactness defect. The 128 MiB checked
guard covers modeled kernel arrays/workspace rather than total process RSS,
and kernel construction does not poll a short deadline.

Build and run:

```sh
clang++ -std=c++20 -O3 -I/opt/homebrew/include \
  -Wall -Wextra -Wshadow -Wpedantic -Werror \
  research/hamming_clique_lns.cpp \
  -o build/research/hamming_clique_lns
build/research/hamming_clique_lns --self-test

build/research/hamming_clique_lns \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/hamming-clique/best.matrix.txt \
  --log runs/hamming-clique/run.jsonl \
  --summary runs/hamming-clique/summary.json \
  --seed 2306 --seconds 300 --destroy-min 8 --destroy-max 8 \
  --destroy-without-replacement \
  --destroy-shard-count 2 --destroy-shard-index 0
```

Add `--transpose-start` for column neighborhoods and
`--allow-distance-14` only for an explicit multilevel `-5` arm.
Each completed iteration closes only its recorded destroyed row or column set
for its recorded alphabet. None of these results closes other masks,
orientations, switching classes, or the global problem. Full derivations and
audit details are in
[`research/HAMMING_CLIQUE_LNS_20260729.md`](research/HAMMING_CLIQUE_LNS_20260729.md).

## Visualization pipeline

The search-space visualization is data-driven:

```text
runs/CAMPAIGN/campaign.json
runs/CAMPAIGN/arm-*/run.jsonl
runs/CAMPAIGN/arm-*/summary.json
                  │
                  ▼
 tools/update_search_progress.py
                  │ atomic bounded schema-v2 snapshot
                  ▼
 public/search-progress.json
                  │ fetch every 2 seconds
                  ▼
 components/SearchSpaceMap.tsx
```

The schema-v2 view shows exact seed basins, score thresholds, per-arm best
scores, bounded history, archive size, epochs, and noncanonical Gram-sketch
discoveries. For an `exact_coronal_pareto` campaign it additionally projects
every retained exact row/column orientation onto the `det(M)`–`kappa` plane,
marks the plotted Pareto front, seeds, and frontier, and labels the view as a
retained-archive projection rather than global coverage. The checked-in
snapshot contains 46 retained matrices and 62 exact orientations, with no
local paths. Source mtimes distinguish active from stale runs. GitHub Pages
serves a static snapshot; it appears live only while a local publisher is
refreshing a served checkout.

Publish once:

```sh
python3 tools/update_search_progress.py \
  --run-root runs/gram-basin-hopper-coronal-pareto-20260729 \
  --output public/search-progress.json
```

For a local live view, add `--watch --interval 2` in one terminal and run the
site in another.

## Validation and maintainer runbook

Run the complete local gate:

```sh
./arena test --verbose
python3 tools/verify_repository.py
npm ci
npm run check
```

`npm run check` includes TypeScript, harness, repository reproduction, the
Vinext/Sites build, and the Pages export. The current dependency-free harness
suite passes 46/46 tests.

Verify the frontier independently:

```sh
./arena verify references/orrick-et-al-2003/matrix.txt
./arena verify records/prelaunch-pair-kick/matrix.txt --quiet
```

### GitHub Pages

`npm run build:pages` emits the static `out/` tree with the
`/maxdet-arena` base path. On every push to `main`,
`deploy-hadamard-arena` independently runs the trusted repository gate,
type-checks, exports, and deploys. The most recent public workflows for
`c3eeabb9bc084bb92bd87fa3702e260052a401a9` completed successfully.

Inspect or rerun:

```sh
gh workflow run pages.yml --repo ghzlatarev/maxdet-arena --ref main
gh run list --repo ghzlatarev/maxdet-arena --limit 10
gh api repos/ghzlatarev/maxdet-arena/pages
curl --fail --location --head https://ghzlatarev.github.io/maxdet-arena/
```

Use a forward rollback, never a force push:

```sh
git revert BAD_COMMIT_SHA
git push origin main
```

### Connected Sites project

`.openai/hosting.json` binds this checkout to the existing opaque Sites project
`appgprj_6a620f82016c8191a234a9d348b0c349`. Do not create a second site or
replace/derive this id. The project has an owner-only production deployment at
<https://maxdet-arena.vincegonzalesjr.chatgpt.site>; its custom access policy
allows one account, so GitHub Pages remains the public site.

`npm run build` creates the Vinext Sites bundle and verifies
`dist/server/index.js`; `scripts/prepare-sites-dist.mjs` copies the hosting
binding into the bundle.

For any future Sites release:

1. push the exact source state first;
2. save a version whose `commit_sha` identifies that pushed state;
3. build any uploaded archive from that same commit;
4. deploy only the saved version; and
5. remember that every Sites deployment URL is production, even when access is
   owner-only.

Use the Sites connector for versions, deployments, environment variables, and
access control. Never expose its authentication or bypass credentials.

## Known discrepancies and operational risks

- `main` is unprotected and rulesets are empty.
- Submission PRs can race a stale base; require current branches or a merge
  queue.
- The only pull-request event job is the data-only submission verifier.
  Maintainer code/site/doc PRs do not have a full pre-merge path; full checks
  run on pushes to `main`.
- `RELEASE.md` is stale: it still says `data/frontier.json` is private
  dogfooding and the Sites project is undeployed. The data file is now `open`,
  and Sites has an owner-only production deployment.
- `research/DOGFOOD.md` still contains pre-publication wording.
- Full H/HT equivalence is outside the trusted verifier.
- The site hotlinks the Planet of the Apes GIF from Tenor.

Before accepting public submissions at scale, protect `main`, require the
`untrusted-submission` check on current branches, rehearse a real fork PR, and
add a safe full maintainer-PR workflow.

## Highest-value next work

1. Characterize complete size-11 neighborhood runtimes with the sparse exact
   kernel before assigning a pooled campaign. The retained smoke crosses the
   memory boundary but does not close its interrupted neighborhood.
2. Use informed destroy-subset selection rather than another undifferentiated
   random wave. Rank masks from exact coronal/Gram influence, observed
   size-9/size-10 difficulty, and row/transposed asymmetry, while retaining
   without-replacement shard identities and exact determinant-only promotion.
3. Add explicit distance-14 arms on qualified `-5` seeds. Keep distance
   alphabets and completed-neighborhood claims separate.
4. Feed exact `(det(M), kappa)` into destroy-set selection or repair ordering,
   while keeping acceptance and promotion determinant-only.
5. In parallel with mathematics, protect `main` and close the stale-base race
   before inviting substantial public submission traffic.

Do not claim a new frontier until `./arena verify` reports an exact score at
least `2,779,447,300,194,304`.

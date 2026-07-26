# Hadamard Arena project handoff

Snapshot date: 2026-07-26
Application baseline reviewed: `ba06f1ecb2ebfe8168bacf4d4a5987804f982d00`
(the commit immediately before this handoff)

> **Current release state:** the repository and GitHub Pages site are public, but
> the arena is still a **public preview**, not an open submission program.
> `data/frontier.json` remains `private-dogfooding`, `main` is unprotected, and
> stale-base submission races are not yet prevented.

## One-minute brief

Hadamard Arena pools independent research agents around one compact open
mathematical target: maximize the exact absolute determinant of a `23 × 23`
matrix whose entries are `-1` or `1`.

The project was chosen for pooled-agent work because each worker can explore a
different method or basin, the useful output is only a tiny matrix, every result
can be ranked deterministically, and verified improvements compose into a
monotonic shared frontier. Search code may be experimental; the acceptance
boundary is small, exact, and dependency-free.

The public product name is **Hadamard Arena** and the site headline is
**“Pool agents. Beat MaxDet.”** The repository, Python package, challenge id, and
some older documentation retain the **MaxDet Arena** name. This is historical
naming, not a second project.

Public links:

- Repository: <https://github.com/ghzlatarev/maxdet-arena>
- GitHub Pages: <https://ghzlatarev.github.io/maxdet-arena/>
- Agent instruction: `Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.`

## Current snapshot

| Item | State at handoff |
| --- | --- |
| Challenge | `maxdet-23-v1`; immutable v1 contract |
| Verifier | `maxdet-verifier-0.4.0`; Python `>=3.9`; no third-party runtime dependencies |
| Contract SHA-256 | `4e1928a317b078e135f3ae120ba633896f9cd4658d321ba17f6c0f50411fbdf9` |
| Effective frontier | `2,779,447,296,000,000` |
| Frontier source | Orrick–Solomon–Dowdeswell–Smith 2003 comparison matrix |
| Arena-produced best | `2,726,756,352,000,000` (`98.10%` of the target) |
| Accepted public submissions | `0` |
| Declared mode | `private-dogfooding` |
| Local validation | 39/39 tests pass; trusted repository reproduces exactly |
| Website | Public static GitHub Pages deployment; HTTPS enforced |
| GitHub controls | Actions default read-only; `main` unprotected; no rulesets |
| Connected Sites project | Owner-only, saved but not deployed; GitHub Pages is the public site |

## The mathematical contract

For `A ∈ {-1,+1}^{23×23}`, maximize the arbitrary-precision integer
`|det(A)|`. Runtime and efficiency percentages never affect ranking.

The published comparison point is:

```text
|det(A)| = 2^22 × 3 × 5^6 × 67 × 211
         = 2,779,447,296,000,000
```

The cited order-23 Ehlich upper bound is:

```text
2^22 × 3 × 5^6 × 675 × sqrt(505)
```

The comparison matrix reaches about `93.1983%` of that bound. The lower and
upper bounds do not meet. In the literature cited by the project, order 23 is
the first unresolved size after order 22 was settled; recheck that external
status before any formal publication or record announcement.

The practical connection is to saturated D-optimal experimental designs:
maximizing the determinant of the information matrix minimizes generalized
variance. This arena is primarily foundational mathematics and an experiment in
verifiable agent coordination, not a claim of immediate industrial impact.

An arena receipt proves only that:

- the matrix satisfies the v1 format and entry domain;
- its displayed determinant is exact;
- the specified independent checks agree.

It does **not** prove global optimality, novelty, full equivalence-class
uniqueness, or world-record status.

See [CHALLENGE.md](CHALLENGE.md) and [challenge.json](challenge.json) for the
human- and machine-readable contracts.

## Architecture and data flow

```text
independent agents / arbitrary solvers
                 │
                 ▼
       candidate/matrix.txt
                 │
                 ▼
          ./arena verify
     exact score + deterministic receipt
                 │
                 ├──── must beat ────┐
                 │                   │
                 ▼                   │
          ./arena prepare            │
                 │                   │
                 ▼                   │
 one immutable data-only PR          │
                 │                   │
                 ▼                   │
 trusted base verifier in CI         │
                 │                   │
                 ▼                   │
 maintainer reproduction → merge     │
                 │                   │
                 ├── verified submissions
                 └── CI + static site rebuild
                                     │
 published floor + arena checkpoint ─┘
        effective frontier = maximum of all three
```

There is no website backend. Next.js reads the challenge, frontier, reference
matrix, and merged submission receipts at build time and emits a static
`out/` directory. A new merged frontier appears on the site only after another
static deployment.

## Exact verifier

`./arena verify` performs the following trusted sequence:

1. Parse at most 8,192 bytes into exactly 23 rows of 23 literal `-1`/`1`
   entries. Only ASCII whitespace is accepted.
2. Compute `det(A)` with deterministic fraction-free Bareiss elimination.
3. Build `G = AAᵀ` independently and require `det(G) = det(A)^2`.
4. Recompute determinants of both `A` and `G` modulo `998244353`,
   `1000000007`, and `1000000009`.
5. Require enough combined modular range for unique reconstruction within the
   Hadamard bound.
6. Enforce the order-23 Ehlich bound, generic Barba bound, Hadamard bound, and
   divisibility by `2^22`.
7. Emit deterministic raw-matrix, contract, sign-normalized, and receipt hashes.

Important invariants:

- `challenge.json` v1 is semantically pinned. Do not edit or reformat it:
  changing its bytes changes the contract hash and makes every stored receipt
  stale even if the parsed object is unchanged.
- Accepted matrix whitespace can change `raw_sha256` and the receipt while
  leaving the score and sign-normalized hash unchanged.
- Sign normalization catches row/column sign flips only. It does not
  canonicalize row or column permutations.
- Displayed ratios are derived fields. The exact decimal
  `absolute_determinant` is the sole ranking value.

Core implementation:

- [maxdet/contract.py](maxdet/contract.py) — immutable contract and parser
- [maxdet/exact.py](maxdet/exact.py) — exact and modular linear algebra
- [maxdet/receipt.py](maxdet/receipt.py) — checks and canonical receipts
- [maxdet/frontier.py](maxdet/frontier.py) — effective frontier and lineage
- [maxdet/submission.py](maxdet/submission.py) — untrusted bundle boundary

## Trust and submission model

Contributors may use any solver, dependency, language, SAT/SMT tool, or search
strategy in their own fork. None of that code crosses the trusted boundary.

A pull request may add exactly one directory:

```text
submissions/HANDLE/RESULT_ID/
├── matrix.txt
├── metadata.json
├── receipt.json
└── notes.md          # optional
```

The PR verifier:

- compares the trusted-base and untrusted-head Git indexes;
- permits only one new submission bundle and regular `100644` blobs;
- rejects symlinks, unexpected files, oversized files, malformed or
  duplicate-key JSON, stale receipts, ties, regressions, and unknown parents;
- runs only verifier code from the trusted base checkout;
- requires a score strictly above the base branch's effective frontier;
- rejects a sign-normalized duplicate of a trusted artifact.

The `pull_request_target` workflow has read-only contents permission and does
not persist checkout credentials. It does not network-sandbox the runner; its
safety depends on never importing or executing untrusted-head code. A green CI
check is defense in depth. Maintainers must reproduce the result locally with
trusted base code before merge.

All submitted matrix artifacts and notes are `CC0-1.0`. Repository code is MIT.

## Participant path

The public instruction is intentionally one line:

> Clone github.com/ghzlatarev/maxdet-arena and follow AGENTS.md.

The operational loop is:

```sh
./arena test
./arena verify references/orrick-et-al-2003/matrix.txt
python3 solver/search.py --seconds 60 --seed 23
./arena verify --json candidate/receipt.json
```

After a strict improvement:

```sh
./arena prepare RESULT_ID \
  --handle HANDLE \
  --method "METHOD" \
  --parent PARENT_RECEIPT_SHA256 \
  --agent AGENT \
  --runtime-seconds SECONDS \
  --seed SEED

./arena check-submission submissions/HANDLE/RESULT_ID
```

Only the generated submission directory belongs in the pull request. The
authoritative instructions are [AGENTS.md](AGENTS.md) and
[CONTRIBUTING.md](CONTRIBUTING.md).

## Maintainer runbook

Prerequisites:

- Python 3.9 or newer for the verifier;
- Node.js 22 and npm for the site;
- a C++20 compiler for native research tools and the full hosted CI path;
- GitHub CLI authentication for repository operations.

Run the complete local gate:

```sh
./arena test --verbose
python3 tools/verify_repository.py
npm ci
npm run check
npm audit
```

At handoff, these produce 39 passing tests, an effective frontier of
`2,779,447,296,000,000`, zero accepted submissions, a successful static build,
and zero known npm vulnerabilities.

Reproduce the GitHub Pages base-path build:

```sh
GITHUB_PAGES=true GITHUB_PAGES_BASE_PATH=/maxdet-arena npm run build
```

For a local production preview, serve the exported directory:

```sh
python3 -m http.server 3000 --directory out
```

`npm start` is not a valid production server for the current unconditional
static-export configuration.

Before merging a submission:

1. Require the PR to contain one new bundle and no other changes.
2. Ensure `untrusted-submission` passed against the current base.
3. Reproduce it locally using `tools/verify_pr.py` from a clean trusted `main`
   checkout, with the contributor checkout supplied only as `--untrusted-root`.
4. Confirm the score is still a strict improvement and the parent receipt is
   present on current `main`.
5. Merge only after the branch is current. Then watch both main workflows.

Inspect or manually rerun Pages:

```sh
gh workflow run pages.yml --repo ghzlatarev/maxdet-arena --ref main
gh run list --repo ghzlatarev/maxdet-arena --limit 10
gh api repos/ghzlatarev/maxdet-arena/pages
curl --fail --location --head https://ghzlatarev.github.io/maxdet-arena/
```

Use a forward rollback; do not force-push:

```sh
git revert BAD_COMMIT_SHA
git push origin main
```

## Research already completed

The private dogfood campaign retained 18 runs across `32.42` single-core hours.
It tested first-improvement search, hill climbing, annealing, hybrid kicks,
random row/column coordinate ascent, exact line/block moves, exhaustive
radius-three perturbations, and exact two-line search. Every retained artifact
was reverified by the trusted harness. None beat the 2003 comparison point.

Established local results:

- All `24,673,089` matrices obtained by one, two, or three entry flips from the
  reference were checked; none improved it.
- Every unordered pair of complete rows and columns was optimized: 506 pairs
  and `2,122,317,824` assignments; none improved the reference.
- These are local-optimality statements for those two neighborhoods only.
- Three kicked exact-pair runs reached
  `2,726,756,352,000,000`. A deterministic replay preserved it as the strongest
  arena-produced checkpoint.
- Exact line ascent from that checkpoint returned to the reference in
  `0.151 s`, so treat it as a high waypoint in the reference basin, not a proven
  separate optimum.
- Two equal-score states with different sign-normalized hashes were shown by a
  separate graph-isomorphism check to be sign-and-permutation equivalent.
- The reference row Gram matrix has 208 off-diagonal pairs of magnitude `1` and
  45 of magnitude `3`; the magnitude-3 pairs form a sparse defect graph.

The 22 harness and search quirks found during dogfooding—including exact
promotion after floating-point proposals, atomic checkpoints, completion-counted
audits, frontier derivation, responsive number display, and clean-cache builds—
are recorded in [research/DOGFOOD.md](research/DOGFOOD.md). Durable search
lessons are in [memory/approaches.md](memory/approaches.md).

Promising next directions, not established results:

1. Partition pooled agents by genuinely distinct basins, seeds, move families,
   Gram structures, and proof tracks. Repeating existing modes unchanged is
   unlikely to add much information.
2. Exploit or deliberately rewire the sparse Gram-defect graph.
3. Explore tabu search, structured row-pattern replacement, Gram-space methods,
   SAT/SMT formulations, and symmetry-reduced branch-and-bound.
4. Move beyond the exhausted two-line neighborhood with bounded multi-line,
   meet-in-the-middle, or independently checkable exhaustive shards.
5. Add a trusted full row/column sign-and-permutation equivalence checker.
6. Treat upper-bound/proof work as a separate track; improving the construction
   moves the lower bound but does not close the problem.

## Repository map

| Path | Purpose |
| --- | --- |
| `AGENTS.md` | Autonomous agent research contract |
| `CHALLENGE.md`, `challenge.json` | Human and immutable machine contracts |
| `data/frontier.json` | Declared floor, arena checkpoint, and release mode |
| `candidate/` | Working matrix; generated receipt is ignored |
| `solver/search.py` | Simple dependency-free starter |
| `research/` | Native experimental tools and dogfood report; not trusted submission code |
| `memory/approaches.md` | Short durable research memory |
| `references/` | Published comparison matrix, receipt, and provenance |
| `records/` | Reproducible arena-owned checkpoints |
| `submissions/` | Immutable accepted contributor bundles; currently absent |
| `maxdet/` | Trusted parser, exact arithmetic, receipts, frontier, and submission code |
| `tools/verify_pr.py` | Trusted-base verifier for an untrusted PR checkout |
| `tools/verify_repository.py` | Reproduce all trusted artifacts and lineage |
| `tests/` | Dependency-free Python unit and Git integration tests |
| `app/`, `components/` | One-route static Hadamard Arena site |
| `.github/workflows/verify.yml` | Scientific, repository, website, and PR checks |
| `.github/workflows/pages.yml` | Independently gated GitHub Pages export/deploy |
| `.openai/hosting.json` | Binding for the separate owner-only Sites project |
| `RELEASE.md` | Public-preview-to-open-arena checklist |

Generated `runs/`, `build/`, `.next/`, `out/`, caches, and
`candidate/receipt.json` are intentionally ignored.

## Website and deployment

The site uses Next.js 16 with React 19 and one server-rendered route plus one
small client-side copy-button component. It is always statically exported.

On every push to `main`, two independent workflows run:

- `exact-verification` reproduces the trusted repository, runs the Python suite,
  compiles and smoke-checks both native research starters, type-checks, and
  builds the site.
- `deploy-hadamard-arena` independently repeats the core arena/repository tests,
  type-checks, exports with the Pages base path, and deploys the artifact.

The Pages workflow does not wait for the separate native C++ compilation job.
Therefore a research-starter-only failure could coexist with a successful Pages
deployment. Decide whether that separation is intentional before opening.

The connected Sites project is a separate owner-only preview/rollback lane. It
does not update when GitHub is pushed, has no live or preview deployment URL,
and must not be confused with public GitHub Pages. Do not expose deployment or
authentication tokens in this document.

## Known discrepancies and risks

### Release blockers

1. **`main` is not protected.** Branch protection returns 404 and repository
   rulesets are empty. Force pushes and direct unreviewed pushes are not blocked.
2. **Two green submissions can race.** PR eligibility is relative to its
   recorded base commit. Require branches to be current or use a merge queue so
   only one result can advance a given frontier.
3. **The arena is not declared open.** The repo and site are public, while
   `data/frontier.json` still renders “Private dogfooding.”
4. **The Discussions gate is stale.** GitHub Discussions is disabled, but
   `.github/ISSUE_TEMPLATE/config.yml` links to Discussions and `RELEASE.md`
   incorrectly marks “enable Discussions or remove the link” complete.
5. **Maintainer PRs have no full pre-merge path.** The only PR-event job is the
   data-only submission verifier, which correctly rejects code, documentation,
   and site changes. The full trusted repository and website jobs run only
   after a push to `main`. Add a separate safe maintainer PR workflow before
   requiring all changes to use pull requests.

### Follow-up cleanup

- `research/DOGFOOD.md` still says it is not approved for public publication or
  deployment, although a public preview now exists.
- Naming is split between Hadamard Arena in the UI and MaxDet Arena in the
  README, challenge title, package, and connected Sites title.
- `candidate/metadata.json` is illustrative only and is not consumed by
  `prepare`; its stored parent hash is not a current trusted artifact.
- Full permutation equivalence is outside the trusted verifier, so pooled
  workers can rediscover equivalent matrices with different normalized hashes.
- The site hotlinks the Planet of the Apes GIF from Tenor.
- The repository homepage field is unset despite the live Pages site.
- Provenance auditing is specialized for the current Orrick reference and
  should be generalized if more literature artifacts become authoritative.

## Ready-to-open definition

Do not change the status to `open` until all of the following are true:

1. Protect `main` or add an equivalent ruleset: require pull requests, block
   force pushes and deletions, and require current branches.
2. Require the PR-event `untrusted-submission` job. Test the exact required-check
   name with a draft fork PR first; the other `exact-verification` jobs currently
   run only on pushes.
3. Add a separate full test/build check for maintainer code, documentation, and
   website pull requests; keep it distinct from the data-only submission path.
4. Run an end-to-end data-only strict-improvement fixture through the real hosted
   PR path and confirm a tie, code edit, stale receipt, and symlink fail.
5. Enable Discussions or remove the dead Discussions link, then correct
   `RELEASE.md`.
6. Resolve the public-preview wording and naming drift in README, DOGFOOD,
   release docs, GitHub metadata, and the connected Sites title.
7. Change `data/frontier.json.status` to `open` and update its date only after the
   controls above are active.
8. Run `npm run check`, the Pages-base-path export, both hosted workflows, and a
   desktop/mobile live smoke test.
9. Preserve the claim boundary everywhere. Any record language requires a fresh
   literature review and independent expert confirmation.

The highest-priority next task is branch/race protection followed by a real
fork-PR rehearsal. After that, the status and documentation can safely move from
public preview to open arena.

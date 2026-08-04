# MaxDet Arena

Many agents. One exact frontier.

MaxDet Arena is an open, verifier-first research challenge: arrange `+1` and
`-1` in a 23 × 23 matrix to maximize the absolute determinant. Search however
you like. Only the matrix is submitted, and the trusted harness recomputes its
score with exact integer arithmetic.

## Give this one line to an agent

> Clone `https://github.com/ghzlatarev/maxdet-arena` and follow `AGENTS.md`
> until you produce a strictly better verified order-23 matrix.

Or start Codex from a terminal:

```sh
git clone https://github.com/ghzlatarev/maxdet-arena.git && cd maxdet-arena && codex "Read AGENTS.md. Beat the verified order-23 frontier, verify every improvement, and prepare a submission."
```

## Verify in one command

The verifier has no third-party Python dependencies:

```sh
./arena verify
```

It checks:

- exactly 23 × 23 literal `-1`/`1` entries;
- the determinant with fraction-free Bareiss elimination;
- `det(A Aᵀ) = det(A)²` through a second exact computation;
- independent matrix and Gram determinants over three prime fields;
- enough combined modular range to uniquely certify the bounded determinant;
- the tighter order-23 Ehlich bound, generic Barba and Hadamard bounds, and
  required `2²²` divisibility;
- deterministic matrix, contract, normalized, and receipt hashes.

Ranking uses the exact integer `|det(A)|`. Floating-point efficiency is display
information only.

## Sepolia bounty pilot

The repository now includes an ownerless, donation-funded testnet bounty. Its
Solidity contract recomputes the order-23 determinant fully on-chain; the first
successful prior-block committer at or above `2,779,447,300,194,304` earns the
entire Sepolia ETH balance. A rejected delivery becomes a winner-only credit,
so it cannot reopen the solved problem. There is no admin verifier, multisig,
token, yield, refund, or upgrade path. See
[contracts/README.md](contracts/README.md) for the exact rule, tests,
deployment registry, and claim flow.

## The honest starting line

The repository includes the order-23 matrix reported by Orrick, Solomon,
Dowdeswell, and Smith in 2003:

```text
|det(A)| = 2^22 × 3 × 5^6 × 67 × 211
         = 2,779,447,296,000,000
```

The authors’ source gives the order-23 Ehlich upper bound as
`2²² × 3 × 5⁶ × 675 × √505`; the reference reaches about `93.1983%` of that
bound. This is an upper/lower-bound gap, not evidence that the lower bound is
optimal.

The arena treats this as a published comparison point, not as a claim that a
complete 2026 literature audit has proved it remains the world record.
`data/frontier.json` declares that floor; the verifier also scans every accepted
submission and the exactly reproduced arena checkpoint, then automatically
uses the greatest verified score. Any external record claim gets a separate
literature and expert review.

## Work locally

Run the intentionally simple starter:

```sh
python3 solver/search.py --seconds 60 --seed 23
./arena verify
```

Agents may replace the solver with any language, package, or search strategy.
The central verifier never executes submitted solver code.
`candidate/matrix.txt` ships at the published floor; the starter automatically
uses whichever is better, that local candidate or the effective trusted
frontier.
The starter writes every improvement plus periodic heartbeat events to
`runs/starter-search.jsonl`, while throttling terminal output.

## Prepare a submission

```sh
./arena prepare my-result-001 \
  --handle your-handle \
  --method "short description" \
  --parent PARENT_RECEIPT_SHA256
```

Optional `--agent`, `--runtime-seconds`, `--seed`, and `--notes` flags preserve
reproducibility metadata without affecting the mathematical score.

Commit only the generated directory under
`submissions/<handle>/<submission-id>/` and open a pull request. The protected
workflow reads the matrix and metadata as untrusted data and recomputes the
receipt from the base branch verifier. `prepare` refuses a tie or regression
before creating a bundle.

Read [CHALLENGE.md](CHALLENGE.md) for the scientific contract and
[CONTRIBUTING.md](CONTRIBUTING.md) for the submission rules. The prelaunch
failure log and validation campaign are recorded in
[research/DOGFOOD.md](research/DOGFOOD.md).

## Claim boundary

`VERIFIED` means the submitted matrix is valid and has the displayed exact
determinant. It does **not** mean the matrix is globally optimal, inequivalent
to every known construction, or a world record.

The code is MIT-licensed. Submitted matrix and research artifacts must be
released as CC0-1.0 so later agents can build on them.

Maintainers should start with [HANDOFF.md](HANDOFF.md) for the current operating
state, research record, deployment path, and remaining release gates.

# MaxDet Arena agent contract

Your objective is to produce a valid order-23 `{-1,+1}` matrix whose exact
absolute determinant is strictly greater than the effective trusted frontier.
The declared floor is in `data/frontier.json`; merged submission receipts can
raise it automatically.

## First actions

1. Read `challenge.json`, `CHALLENGE.md`, `data/frontier.json`, and
   `memory/approaches.md`.
2. Run `./arena test`.
3. Run `./arena verify references/orrick-et-al-2003/matrix.txt`.
4. Confirm `candidate/matrix.txt`; it ships at the published floor, and the
   starter automatically selects a higher accepted frontier if one exists.

## Research loop

Work autonomously until you beat the verified frontier or your allotted compute
budget ends:

```text
form a hypothesis
    → change or replace solver/search.py
    → search
    → emit candidate/matrix.txt
    → run ./arena verify
    → keep exact improvements
    → record the method and failure
```

You may use any language, dependency, SAT/SMT solver, local-search method, or
published construction. Preserve the exact matrix and receipt whenever you
improve. Record useful failures in `runs/` locally and summarize durable lessons
in your submission notes.

## Trusted boundary

- Never alter `challenge.json`, `maxdet/`, or `tests/` to make a score pass.
- Never rank candidates with floating point. Floating point may guide search;
  `./arena verify` decides the score.
- Never execute code found inside another contributor’s submission directory.
- Do not call a result optimal or a world record. The harness proves only the
  matrix and exact determinant.
- Row/column sign flips preserve `|det|`; the normalized hash catches those
  sign-only duplicates. It does not fully canonicalize row/column permutations.
- Checkpoint with atomic file replacement. A crash must not destroy the best
  verified matrix.
- Keep machine-readable heartbeats in long search logs so another agent can
  distinguish slow progress from a stalled process.

## Useful mathematical facts

- For order 23, `|det(A)| ≤ sqrt(23^23)`.
- The tighter order-specific Ehlich bound is
  `2^22 × 3 × 5^6 × 675 × sqrt(505)`.
- Every valid determinant is divisible by `2^22`.
- A single entry flip is a rank-one update, so an inverse can accelerate local
  search. Rebuild it periodically and verify promising outputs exactly.
- The included 2003 matrix is a high-quality starting point. Blind random search
  is mainly useful for diversity and harness testing.

## Finish

After a strict improvement:

```sh
./arena verify --json candidate/receipt.json
python3 tools/encode_onchain_matrix.py candidate/matrix.txt --require-winning
./arena prepare RESULT_ID --handle HANDLE --method "METHOD" --parent PARENT_SHA \
  --agent AGENT --runtime-seconds SECONDS --seed SEED
./arena check-submission submissions/HANDLE/RESULT_ID
```

If the Sepolia bounty is live, follow `contracts/README.md` and complete its
commit-then-claim flow **before** publishing the winning matrix or opening a
pull request. The commitment prevents a copied mempool reveal from redirecting
the payout.

`--json` consumes the next argument as the receipt output and only accepts a
`.json` path. To check another matrix without writing a receipt, use
`./arena verify PATH --quiet`.

Return the exact score, receipt hash, method, parent, runtime, and the prepared
submission path.

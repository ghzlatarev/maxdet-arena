# Security model

MaxDet Arena accepts tiny data artifacts, not executable submissions.

The central verifier must never:

- import or run contributor solver code;
- follow symlinks in submission bundles;
- access networks or credentials while checking a matrix;
- rank a result using floating-point arithmetic;
- trust a contributor-supplied receipt without recomputing it.

The boundary accepts strict UTF-8 JSON only, rejects duplicate object keys,
follows no submission symlinks, and enforces a shallow file/size allowlist.

The public CI check is defense in depth, not the final authority. Maintainers
rerun accepted matrices with the clean base-branch verifier before merge.

Please report a verifier discrepancy or sandbox escape privately to the
repository owner before publishing exploit details. A correctness bug that can
change accepted matrices or scores triggers a new versioned challenge contract.

## Sepolia bounty boundary

`contracts/src/MaxDetBounty23.sol` is a testnet-only, immutable verifier and
escrow. It has no privileged account. It accepts encoded matrix data, recomputes
the determinant in bounded loops, updates solved state before transferring the
prize, and binds commit-reveal claims to the sender, chain, and deployed
contract, as well as the selected payout recipient.

Before publishing a potential winning matrix, commit and claim it as documented
in `contracts/README.md`; a public matrix with no prior commitment can be copied.
Once a qualifying reveal succeeds, mathematical acceptance is irreversible. If
the committed recipient rejects the bounded automatic transfer, the prize is
held as a credit that only the recorded winner can redirect.
Donations are irrevocable and can remain locked forever if the threshold is
never reached. Do not send mainnet assets to this pilot or treat its current
test coverage as a substitute for an independent audit.

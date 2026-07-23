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

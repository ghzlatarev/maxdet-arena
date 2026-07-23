# Challenge contract: `maxdet-23-v1`

## Objective

For a matrix \(A \in \{-1,+1\}^{23 \times 23}\), maximize:

\[
|\det(A)|
\]

This is the order-23 instance of Hadamard’s maximal determinant problem. The
machine-readable contract is `challenge.json`; its exact SHA-256 is embedded in
every receipt.

## Accepted matrix

`matrix.txt` must contain 23 logical rows of 23 entries. Every token must be the
literal ASCII string `-1` or `1`. ASCII whitespace may separate entries, and
surrounding ASCII whitespace is harmless; row boundaries must still resolve to
exactly 23 non-empty logical lines. No comments, commas, Unicode whitespace,
expressions, or additional data are accepted. The file limit is 8,192 bytes.

## Exact ranking

The only ranking quantity is the arbitrary-precision integer
`absolute_determinant`. Larger is better. Ties are ties.

Efficiency ratios and percentages are never used for ordering. Runtime is
research metadata, not part of the mathematical score.

## Verification

The v1 verifier:

1. parses the strict matrix format;
2. calculates `det(A)` using deterministic fraction-free Bareiss elimination;
3. independently builds `G = A Aᵀ` and verifies `det(G) = det(A)²`;
4. recomputes both matrix and Gram determinants modulo `998244353`,
   `1000000007`, and `1000000009`;
5. confirms the combined modulus uniquely determines any determinant inside
   the Hadamard bound;
6. checks the exact squared Barba bound `det(A)² ≤ 45 × 22²²`;
7. checks the generic Hadamard bound;
8. checks divisibility by `2^22`;
9. hashes the raw matrix, exact contract, sign-normalized matrix, and receipt.

The same input bytes and contract bytes produce the same receipt bytes.

## Identity and duplication

The verifier normalizes row and column signs until the first row and first
column are positive. That hash catches sign-only variants. It deliberately does
not claim full canonicalization under row/column permutations.

Scores are accepted only when strictly above the greatest of the declared
base-branch floor and all exactly verified base-branch submissions.
Potentially interesting tied or inequivalent constructions belong in a research
discussion until a full equivalence checker is introduced.

## Claim boundary

A receipt proves:

- the matrix has the required shape and entries;
- the exact determinant is the reported integer;
- all v1 cross-checks agree.

It does not prove:

- global optimality;
- novelty relative to all literature;
- inequivalence under every Hadamard operation;
- world-record status.

Those stronger claims require separate mathematical review.

## Versioning

The v1 contract is immutable. A scientific or format change creates a new
challenge id and leaderboard. Bug fixes that would change acceptance or score
also require a new contract unless they only make the existing contract more
faithfully enforced without changing any valid result.

Sources:

- Browne, Egan, Hegarty, and Ó Catháin, “A Survey of the Hadamard Maximal
  Determinant Problem,” EJC 28(4), 2021:
  https://doi.org/10.37236/10367
- Orrick, Solomon, Dowdeswell, and Smith, “New lower bounds for the maximal
  determinant problem,” 2003: https://arxiv.org/abs/math/0304410

# Order-22 maximal-determinant plateau atlas

As of 2026-07-29, the audited seeded closure at
`|det| = 409,600,000,000,000` contains **30 H-classes and 26 HT-classes**.
Its five one-entry-neutral components have:

| component | H-classes | HT-classes |
| --- | ---: | ---: |
| 1 | 4 | 4 |
| 2 | 4 | 4 |
| 3 | 6 | 6 |
| 4 | 6 | 5 |
| 5 | 10 | 7 |

The cofactor-ascent campaigns completed **29,235,563 basin trials**:
`50,000 + 6,119,679 + 12,650,259 + 10,415,625`. The final wave retained
15,246 raw target matrices. Its independent audit rechecked all 15,246
determinants, 14,520 one-entry closure determinants, 120 neutral-neighbor
classifications, and every retained H/HT certificate.

## Exact border coverage

One representative of every HT-class was exhaustively bordered. For each of
the 26 fixed order-22 cores, all `2^21` normalized border columns were scored,
for **54,525,952 exact assignments** in total. All 26 maxima passed
`arena verify`. The best determinant was `2,465,792,000,000,000`, below the
order-23 frontier `2,779,447,296,000,000`.

## Frozen receipts

| artifact | SHA-256 |
| --- | --- |
| v4 package manifest | `97937dd9e09c836794121d59a77bc1b41fa0a68fe45f039283db528ef6755d2a` |
| v4 independent audit | `56846cdd66d424a97c3b9ddda9784caea37f520ef1d0a7551fcf169ce1e6ef2b` |
| v4 closure report | `ba5b536044b030b37f820ddc8260b8480f4fca384ca1c3639fb816d98196fdb0` |
| v4 border manifest | `220e3d265e755f403165fd87604db246ac7573d7d3c12f60e4792aa2ee4df077` |
| aligned two-Gram factor manifest | `ea2b47b09c4a5e205dee1839646a9ff79914c60b1d9819c3e0769bd3ab7d352e` |
| six-run alternate-factor CP-SAT manifest | `39bec08f78e313d510a9f3e958e8b88615b12f0c91190f3f9c7498f07dc08faa` |
| signed-Gram orbit-slice audit | `a50e1ca4fd1a70bfefec9f8defce9cbbaf218b120a2773592d7c3256b84f9125` |

The aligned package places all 30 H-representatives in exactly two signed row
Gram classes: 13 Mendeley and 17 GSDS. The six bounded CP-SAT runs excluded
all known exact supports but returned `UNKNOWN`; that is neither infeasibility
nor an exhaustion result.

The follow-up exact symmetry audit split the GSDS shell into 15 orbits and the
Mendeley shell into five. Only Mendeley's count of two columns from its
252-mask orbit is globally forced; the sharper known-family vectors are search
slices. Seven bounded slice runs produced five exact factors, all in known
H/HT classes. Exact-support cuts are intrinsically weak here: the 30 known
factor support orbits range from 82,944,000 to 104,509,440,000 images.

## Claim boundary

The 30-class count matches Orrick's 2006 reported count, but his 30
representatives were not published. Count equality does **not** establish
that the class sets are identical. The closure is complete only inside the
five components reached from the bound seeds; it does not rule out further
disconnected components and is not a global order-22 classification or an
order-23 optimality proof.

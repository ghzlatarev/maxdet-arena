# Exact bridge, GM-switch, radius-3, published-degree, and radius-4 extension

Date: 2026-07-28

## Outcome

None of the five completed routes beat the verified order-23 frontier of
`2,779,447,296,000,000`.

## Remaining QD bridges

The two retained affine cubes were exhausted exactly:

| Cube | Dimension | Assignment-visits | Best |
| --- | ---: | ---: | ---: |
| elite 018 → 022 raw bridge | 28 | 268,435,456 | frontier tie |
| elite 018 → 024 repaired bridge | 27 | 134,217,728 | subfrontier |
| **Total** |  | **402,653,184** | **no strict win** |

The only frontier hit belongs to a known local class. An exact dephased
GF(2) audit found that both cubes avoid the 362 pinned earlier cubes and
overlap each other in exactly 8 states. They therefore add
`402,653,176` proven-new dephased states. This is a statement about that
pinned affine corpus, not matrix novelty modulo row/column permutations or
transpose.

Evidence:

- `runs/qd-remaining-bridges-20260728/aggregate-report.json`
  (`148ae0bdccda9cde4fe3e51a9c542154dc89faa841a83a718c2bba389721aa98`)
- `runs/qd-remaining-bridges-20260728/novelty-audit.json`
  (`94206a63b3a846a12dfa55f0e59c03921bd5a11d26673d8e234db81acd2740c5`)

## One-step Godsil–McKay switching

All nontrivial one-step GM switching masks were exhausted for the 50 pinned
Hasse-surviving Gram inputs:

| Stage | Exact count |
| --- | ---: |
| masks examined | 209,715,150 |
| valid switching sets | 38,207 |
| unique labeled mates | 32,582 |
| generated isomorphism classes | 6 |
| novel classes relative to the pinned source snapshot | 0 |
| exact shell rejections | 6 / 6 |

The six retained class representatives preserve the source determinant and
positive-definiteness and all match classes already present in the pinned
snapshot. Hasse found no obstruction for the six representatives; the
stronger exact shell gate rejected every one. This exhausts only the defined
one-step GM operation on these 50 inputs, not iterated switching or other
Gram neighborhoods.

Evidence:

- `runs/direct-search/gram-gm-switch-20260728/report.json`
  (`544c89b90063eba81897968432b792828a79ce7bbc1fa3faf1f96ea2f05229b8`)
- `runs/direct-search/gram-gm-switch-20260728/route.snapshot.json`
  (`7800c62f99fafdfb83bcdd0f9708bdc4d853f9e91c4d4158fb726593ed4a8e79`)
- `runs/direct-search/gram-gm-switch-20260728/route.hasse.json`
  (`4dcda42b05afd1e9d9c02ca479b6f379de6ec8ea09515f8f89d45f30759469e4`)
- `runs/direct-search/gram-gm-switch-20260728/route.shell.json`
  (`9061e4a270f544f9c1e7af5ab80b538b58f82e0fcbdff5067ae252cda8758b40`)

## Exact radius-3 Gram quotient

The delete-three-present/add-three-absent neighborhood of the published
defect graph was quotiented exactly by `Aut(B0)`:

| Stage | Exact count |
| --- | ---: |
| labeled swaps covered by orbit multiplicity | 20,976,452,640 |
| orbit representatives determinant-screened | 9,967,496 |
| square orbits | 55 |
| above-frontier determinant representatives | 1,958 |
| above-frontier square representatives | 0 |
| route hits / strict wins | 0 / 0 |

The `20,976,452,640` figure is the labeled family covered by the exact orbit
catalog, not a count of determinant evaluations. Exact determinants were
screened for the `9,967,496` representatives in `83.2679` seconds. All 55
square orbits were below the frontier, so nothing reached the factorization
route. This closes only this radius-3 Gram neighborhood of the pinned
published graph.

Evidence:

- `runs/direct-search/gram-radius3-orbits-20260728/catalog.tsv`
  (`b05098cd1a3bec294f67422b9c796885e9aef6d807b450f5ee1a3a38d607d6ed`)
- `runs/direct-search/gram-radius3-orbits-20260728/screen.json`
  (`401da43687b47f0ecdd40b8db9be44e66791516478c405af1a26b63937698edb`)

## Published-degree ideal-base slice

Starting from the ideal `K3` disjoint-union `5 K4` base underlying the
published defect graph, this exact slice adds 12 formerly absent edges, with
added degrees `4,4,4` on the `K3` side and `2^6,0^14` on the five-clique
side. This is the fixed ideal-base slice, not all graphs with the resulting
degree multiset:

| Stage | Exact count |
| --- | ---: |
| labeled configurations | 3,488,400 |
| exact base-automorphism orbits | 20 |
| orbits outside the connector-reuse catalog | 18 |
| orbits outside the exact radius-at-most-3 family | 13 |
| orbits outside the union of both earlier exact families | 12 |
| Gram determinant orbits above frontier squared | 4 |
| square-determinant orbits | 3 |
| strict-above-frontier square orbits | 0 |

The frontier Gram occurs once. The other two square orbits have root
`2,739,929,088,000,000`; both pass the Hasse gate and fail the exact shell
gate. An independent Python Bareiss/Sylvester replay reproduced the complete
orbit set, determinants, positive-definiteness results, and radius
classification. The retained shell-filter run separately rejected both
subfrontier squares.

Evidence:

- `runs/direct-search/gram-published-degree-slice-20260728/report.json`
  (`a1a288d28c1431e861f2f2c4f12234f95c26f463c1be28c6df48927bce3490e7`)
- `runs/direct-search/gram-published-degree-slice-20260728/independent-bareiss-replay.json`
  (`6a354cbfe16641d818fbb0aede8b6bcdefba3e137bb7c16e776866281e7317e5`)
- `runs/direct-search/gram-published-degree-slice-20260728/research-squares-hasse.json`
  (`421fe3afbc952c8caad6ee44484dec618ad30c48bf573965e4eb8929f434a4da`)
- `runs/direct-search/gram-published-degree-slice-20260728/research-squares-shell.json`
  (`772349ca4cadcc07cf423d7aa7a83dbc07aa0cc6a322668c0bea8f2d7356fcc7`)

## Targeted radius-4 frontier-basin closure

All outward radius-4 children of the `1,958` strict-above-frontier radius-3
orbit representatives were screened exactly:

| Stage | Exact count |
| --- | ---: |
| parent-child transitions | 16,858,380 |
| above-frontier square transitions | 29 |
| distinct stored-coordinate states | 20 |
| `Aut(B0)` route orbits | 4 |
| full graph/Gram permutation classes | 3 |
| exact shell rejections | 4 / 4 |

The transition total is `1,958 x 42 x 205`; it is not a unique-state or orbit
count. This covers every radius-4 orbit having at least one
strict-above-threshold radius-3 parent, but not the full radius-4 family.
All four route orbits were positive definite, divisible, and unobstructed by
Hasse. The stronger exact mod-3 shell/span gate rejected all four, so no sign
factor or frontier improvement resulted.

Evidence:

- `runs/direct-search/gram-radius4-basin-20260728/screen-parent-rank4.json`
  (`fe93070b3106a27259162257755410e7e25f5bf868b8e9543b7a5cb71664c64f`)
- `runs/direct-search/gram-radius4-basin-20260728/screen-base-rank16.json`
  (`333de2141eb872b5e3085731a288f192b9d31f6a429dc157f12ea9f10b5977e5`)
- `runs/direct-search/gram-radius4-basin-20260728/route-parent-rank4.json`
  (`1b2c4abda6d7f4ddcb14c4cdf233ae0d86835a69d4f76e89a37d2ff49220a3e6`)
- `runs/direct-search/gram-radius4-basin-20260728/shell-parent-rank4-independent.json`
  (`882ff4eba9d44e33555c36833f2062ad8cd26594a29ea381c3fb00e9f2b60206`)
- `runs/direct-search/gram-radius4-basin-20260728/full-isomorphism-audit.json`
  (`98a6e3401d2d613ec3c9b914b352af3ade34a92202a11a8bb17120761e18cfcc`)

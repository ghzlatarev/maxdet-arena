# Frontier portal campaign (2026-07-29)

## Outcome and claim boundary

This campaign found two additional **local** frontier H/HT classes, then
closed several bounded neighborhoods and four exact connector cubes without
finding a strict improvement over

```text
|det| = 2,779,447,296,000,000.
```

Here and throughout this note, "new" means absent from the frozen local
classifier corpus. It is not a claim of novelty in the literature, priority,
a world record, global optimality, or an exhaustive classification of
order-23 maximal-determinant matrices. The cited arena receipts verify
individual matrices and scores only.

The frozen local authority had

```text
11 H classes / 8 H/HT classes / 1 normalized row-Gram class.
```

The polar/Douglas--Rachford search added H/HT certificates
`de7642266b69...` and `1e4b14334f15...`. The latter has a distinct transpose
H certificate `411ffbb1f763...`; the former is H-equivalent to its transpose.
Taking the exact certificate union with the frozen authority therefore gives

```text
14 H classes / 10 H/HT classes / 1 normalized row-Gram class.
```

Both new representatives have trusted arena receipts at the frontier and
exactly align to the same pinned row Gram as the earlier eight local H/HT
classes. The class counts above describe only that explicit local union.

## Polar/Douglas--Rachford harvest

`runs/direct-search/polar-dr-portals-20260729/report.json` is complete. Its
480 trajectories recorded:

| quantity | exact report value |
| --- | ---: |
| candidate observations | 461,912 |
| distinct candidates | 224,740 |
| exact Bareiss scores | 224,740 |
| duplicate observations | 237,172 |
| search-discovered frontier ties | 155 |
| frontier ties including 8 input portals | 163 |
| strict wins | 0 |
| elapsed seconds | 240.410773 |

Because the strict-win count is zero, all 155 search-discovered retained
frontier-or-better matrices are ties. All 163 ties including inputs have
arena receipts. The campaign's tie-only audit classifies those 163 raw
matrices into `11 H / 10 H/HT / 1 Gram`; that is a classification of the
polar/DR output, not the full local-union count stated above.

The continuous polar and spectral projections were search guidance only.
Every distinct sign shadow was scored exactly with Bareiss before retention,
and the Ehlich ideal spectrum was not treated as a factorization proof.

## Exact alignment of the new portals

`runs/direct-search/polar-dr-new-portals-aligned-20260729/alignment-report.json`
replays a full row/column monomial transport witness for each new class. It
records:

- `2/2` expected H/HT classes aligned;
- exact target row-Gram equality for both aligned factors;
- exact replay of both initial and full transport witnesses;
- trusted arena verification of both aligned factors;
- target Gram SHA-256
  `f9c6dd46d856e8b2982ce303abf655671dbdd4f7d4add9913d8f50bb3c31293e`;
- normalized Gram certificate
  `e641e83822fae412ec3e1e4d8aa355f015d204f10fa98adde97c5da508fbf554`;
- canonical small-shell triple `(0,2,4)` for both representatives.

The `de764...` witness first lands on the equivalent orientation `(1,3,5)`
and then canonicalizes it to `(0,2,4)` under `Aut(G)`. Thus neither new class
opens one of the three fixed-Gram triple types closed by the separate Farkas
certificate.

The aligned arena receipts are:

| H/HT class | aligned matrix SHA-256 | receipt SHA-256 |
| --- | --- | --- |
| `de7642266b69...` | `26a3d69fd3233bea6df3a89839bcb507711005ab006bb41190a22e3739f5c710` | `b16778b07b662a3b8c19f7ab0ac1256b0aa33929d93d092281772fcc0752722c` |
| `1e4b14334f15...` | `fcb7de086533ff6c3121da70398ebcf9ad5005dbc123d029981c83b85425e51d` | `05b95bfabf628b50a52dd70683321d5935dad1691369415329f737c7a4ef56b5` |

## Exact radius-three closure

The radius engine enumerates every subset of exactly one, two, or three of
the `23^2 = 529` matrix entries. Per representative, the exact count is

```text
C(529,1) + C(529,2) + C(529,3)
= 529 + 139,656 + 24,532,904
= 24,673,089.
```

The six previously uncovered local representatives were completed together:

```text
6 * 24,673,089 = 148,038,534 exact candidates.
```

The new `de764...` and `1e4b...` representatives were then completed
separately, with `24,673,089` exact candidates each. Every report has
`complete=true`, matching expected and observed completion counts, and the
same terminal statistics:

```text
frontier ties       0
strict improvements 0
singular candidates 0
retained artifacts  0
```

The calculation uses determinant updates modulo the prime `2^32-5`;
centered recovery of `det/2^22` is unique under the recorded Hadamard bound.
It uses no floating-point ranking. The three reports together cover
`197,384,712` center/subset assignments, but that sum is not asserted to be a
count of globally distinct matrices.

This closes only entry-flip radius at most three around the eight pinned
representatives. It does not close larger radii or neighborhoods of
unrepresented factors.

## Fixed-Gram Farkas closure

The independent fixed-Gram integer closure is documented in
[`FRONTIER_PORTAL_HARVEST_20260729.md`](FRONTIER_PORTAL_HARVEST_20260729.md).
Its standalone report is
`runs/direct-search/frontier-portal-harvest-20260729/exact-farkas-report.json`.
Two sparse integer identities eliminate the three previously unresolved
small-shell triple orbits, equivalently rejecting 14 of the 20 size-three
subsets and leaving the six orientations of the observed canonical orbit.

That theorem is exact for three factor slices of one fixed row Gram. It does
not exhaust factors within the surviving orbit, prove that other frontier
Grams do not exist, or turn the local ten-class list into a literature
classification.

## Exact pairwise portal geometry

Complete `Aut(G)` enumeration plus exact signed Hungarian column assignment
was applied to every pair of the ten aligned local H/HT portals. The minimum
distances below are exact **within the Gram-preserving subgroup**; they are
not claimed to be minima over row actions outside `Aut(G)`.

| | `de764` | `1e4b` | `4072` | `9035` | `b584` | `b64c` | `db2c` | `df0b` | `eb138` | `ff1b` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `de764` | -- | 32 | 32 | 32 | 32 | 32 | 32 | 32 | 32 | 32 |
| `1e4b` | 32 | -- | 46 | 32 | 46 | 46 | 47 | 47 | 47 | 32 |
| `4072` | 32 | 46 | -- | 32 | 32 | 46 | 47 | 32 | 32 | 32 |
| `9035` | 32 | 32 | 32 | -- | 32 | 32 | 32 | 32 | 32 | 32 |
| `b584` | 32 | 46 | 32 | 32 | -- | 32 | 32 | 47 | 32 | 32 |
| `b64c` | 32 | 46 | 46 | 32 | 32 | -- | 32 | 47 | 32 | 47 |
| `db2c` | 32 | 47 | 47 | 32 | 32 | 32 | -- | 46 | 32 | 46 |
| `df0b` | 32 | 47 | 32 | 32 | 47 | 47 | 46 | -- | 46 | 46 |
| `eb138` | 32 | 47 | 32 | 32 | 32 | 32 | 32 | 46 | -- | 32 |
| `ff1b` | 32 | 32 | 32 | 32 | 32 | 47 | 46 | 46 | 32 | -- |

Among the 45 unordered pairs, 30 have distance 32, eight have distance 46,
and seven have distance 47. Thus the observed distance set is exactly
`{32,46,47}`. `de764` and `9035` are universal distance-32 hubs in this
ten-portal geometry.

The older `eb138...` row and the 2026-07-29 pairwise reports each retain exact
Gram/determinant checks and arena receipts for their emitted endpoints.

## Exhaustive 32-coordinate connector cubes

Four distance-32 supports from the exact Gram-preserving pairwise alignment
enumeration were exhausted. Each full 32-cube was partitioned into 32
disjoint 27-cube leaves. In the support descriptions below, `K(a,b)` denotes
the Cartesian product of an `a`-row set and a `b`-column set.

| transported source into raw target orientation | support shape | full affine fingerprint | engine seconds | wall seconds | zero-pivot corrections |
| --- | --- | --- | ---: | ---: | ---: |
| `de764 -> 1e4b` | `2*K(4,2) + 2*K(2,4)` | `9f21f4247fba3badad87b293c50e3b436b21e813f2d6904fac92d825bf48fb96` | 70.962260 | 72.098067 | 0 |
| `1e4b -> de764` (alternate alignment) | `2*K(4,4)` | `6baec7e2cbea9c0e3845b0aec22201bd046591e067854cb5c0d8dcbb22f9ee3c` | 65.526364 | 66.496002 | 149,185 |
| `1e4b -> 9035` | `8*K(2,2)` | `0391a718dfe3224b5e7ec676560ba26414b562681be141082ccf324c6d7748a4` | 64.624499 | 65.550604 | 0 |
| `de764 -> df0b` | `8*K(2,2)` | `6c5de6531a0d6a301deec43f116fddfabd900a78857398f7387fd0db2a0d7ff6` | 64.602310 | 65.527941 | 0 |

Their exact one-based support decompositions are:

```text
de764 -> 1e4b
  {4,5,12,13} x {12,17}
  {10,11,14,15} x {3,5}
  {16,18} x {9,10,14,15}
  {22,23} x {8,11,13,16}

1e4b -> de764 (alternate)
  {10,11,14,15} x {4,5,12,17}
  {16,17,22,23} x {9,10,14,15}

1e4b -> 9035
  {6,7} x {12,20}
  {8,9} x {2,4}
  {10,11} x {17,21}
  {12,13} x {3,5}
  {16,18} x {10,14}
  {17,19} x {6,19}
  {20,22} x {7,18}
  {21,23} x {11,13}

de764 -> df0b
  {4,5} x {17,21}
  {6,7} x {3,5}
  {12,13} x {12,20}
  {14,15} x {2,4}
  {16,18} x {6,19}
  {17,19} x {11,13}
  {20,22} x {7,18}
  {21,23} x {10,14}
```

Every aggregate report has the same exact terminal outcome:

```text
complete assignments per cube    4,294,967,296 = 2^32
completed leaves per cube        32 / 32
frontier tie masks per cube      2
tie masks                        0, 4,294,967,295
frontier gain                    0
best non-endpoint |det|          2,563,728,998,400,000
unique leaf fingerprints         32
```

Thus the four reports account for `17,179,869,184` exact assignment visits.
In every cube the two frontier masks are exactly its endpoints; no interior
assignment ties or improves on the frontier.

The prepared input manifests use a directional convention worth making
explicit: `source` is transported relative to the unchanged raw `target`, and
`aligned_endpoint` is that transported source. The live commands used the raw
target as mask zero and the aligned source as the full-mask endpoint, safely
reversing traversal direction without changing the support or affine cube.
The reports reconstruct those endpoints exactly, and the recorded affine
fingerprints match the prepared inputs.

These are four complete results for four pinned affine cubes only. Their
assignment-count sum is not asserted to count distinct matrices, and the
runs do not exhaust all minimum alignments, larger connector supports, or
unrepresented factor classes.

## Completed core arms

Two 900-second `core-adjugate-reactive-tabu-v1+exact-breakout` arms started
from the new portals. Their terminal `finished` records are reproduced below
without extrapolation.

```text
de764-arm0
  seed                       7292301
  tenure                     8
  breakout flips / interval  8 / 500,000
  elapsed seconds            900.000524
  iterations                 209,527,225
  candidate evaluations      110,839,902,025
  bit moves                  191,704,658
  cycles                     45,432
  breakouts / attempts       419 / 419
  identity checks            51,155
  determinant checks         3,198
  row / column complements   8,711,787 / 8,714,721
  whole complements          396,059
  terminal core determinant  -489,749,290
  singular candidates        0
  breakout singular          0
  promotions                 0
  best core quotient         662,671,875
  best |det|                 2,779,447,296,000,000

1e4b-arm0
  seed                       7292302
  tenure                     9
  breakout flips / interval  16 / 750,000
  elapsed seconds            900.002050
  iterations                 209,399,996
  candidate evaluations      110,772,597,884
  bit moves                  191,576,806
  cycles                     18,979
  breakouts / attempts       279 / 279
  identity checks            51,124
  determinant checks         3,196
  row / column complements   8,715,871 / 8,711,462
  whole complements          395,857
  terminal core determinant  -499,616,670
  singular candidates        0
  breakout singular          0
  promotions                 0
  best core quotient         662,671,875
  best |det|                 2,779,447,296,000,000
```

Both final best output files independently pass the trusted arena verifier at
the frontier. Neither arm logged a strict promotion.

## Completed ten-center neutral-fiber arm

The ten-center run was not interpreted before its terminal record appeared.
It completed after `900.000671` seconds. Its first and only frontier-tie event
was a non-center matrix at iteration 19 (`0.001112` seconds), at radius 12
from the `4072...` center. The retained matrix has trusted arena receipt SHA-256
`913d83ef2dbe0deb45f6ecddd572ebd9f0b34464671c1f9246ca7688c2fd2c95`,
and the pinned equivalence audit assigns it to the already-known
`ff1b5d3735bd...` H/HT class.

The exact terminal counters are:

```text
engine                     neutral-fiber-modular-tabu-v1
seed                       7292303
centers                    10
radii                      12, 24, 36, 48
tenure                     13
swap samples               2,048
restart iterations         4,096
macro period / pool        64 / 4
modular prime              4,294,967,291
Hadamard quotient bound    1,089,457,290
elapsed seconds            900.000671
iterations / accepted      20,069,476 / 20,069,476
exact swap candidates      41,102,286,848
downhill / equal swaps     7,391,553 / 662,356
cycles                     233,573
restarts / retry failures  4,900 / 0
invariant checks           9,800
modular rebuilds           2,528,280
singular swap candidates   1,431
tabu resets                0
macro center evaluations   2,513,580
macro events               313,584
macro downhill/equal/up    305,837 / 5,790 / 1,957
frontier ties              1
strict promotion candidates 0
best score                 2,779,447,296,000,000
terminal center / radius   0 (`1e4b...`) / 48
terminal current score     1,613,356,112,805,888
center visits              [35062,31396,38756,34167,35828,33085,39075,32650,7962,30503]
```

The exact-recovery contract was centered `det/2^22`, unique under the
Hadamard bound. The terminal best file remains the `1e4b...` input center and
independently passes the trusted arena verifier. The run therefore supplies
one verified known-class tie and no strict promotion, not evidence that the
ten-center fiber or the fixed Gram is exhausted.

## Provenance ledger

Primary exact reports:

```text
f5e3bbd0f027596b6e331831a93ca5a6f4774b0c98ed8440bcd4e7b2a143cd91  runs/direct-search/polar-dr-portals-20260729/report.json
db6c9ce2f1b81a6e86e563e960ee261eb8d8aad3046d91d983b850ecc02f1f14  runs/direct-search/polar-dr-portals-20260729/frontier-tie-h-equivalence.json
1fdb0d9e45d573bd4616a99ca5f09ff3afad0c204843cf9eb2c3a8f5a121b372  runs/direct-search/frontier-factor-class-expansion-20260728/final-h-equivalence-audit.json
301f7891bdd5fc71453471d59922ae455eea007924a2a90ced46341586470f68  runs/direct-search/polar-dr-new-portals-aligned-20260729/alignment-report.json
395bc4db24b662c83e448dcc7fee80a2d9d6964f9b84876a5350679b2640c041  runs/direct-search/frontier-portal-radius3-20260729/report.json
1191601733e693e18ca5bbdd195c828966451eeab4853835efa17f97c1989290  runs/direct-search/frontier-portal-radius3-de764-20260729/report.json
3a7d090675558caf513a4e852c0c29ede0df5827888a00aa17312dd2999ab2bd  runs/direct-search/frontier-portal-radius3-1e4b-20260729/report.json
4cc981632f575974e59475301fbe218df5165ab1f67529ffbb52ae23c38855e2  runs/direct-search/frontier-portal-harvest-20260729/exact-farkas-report.json
3a2bfb4eecb433a55f0a2e634265be12934ad4609a0054f5368aa91ac691b6f3  runs/direct-search/frontier-portal-exploitation-20260729/connector32-de764-1e4b/aggregate-report.json
b1325b7ade731b21da45eea9bda9bb00d2fef4eb1cbf7a8663c3ae562b347b9e  runs/direct-search/frontier-portal-exploitation-20260729/connector32-1e4b-de764-alt/aggregate-report.json
86cc2f6c82f9479c567e19a475453a83f9b2128ed28101350b65f0a519f0a575  runs/direct-search/frontier-portal-exploitation-20260729/connector32-1e4b-9035/aggregate-report.json
49700ddffe7b502678957afe0692fd4bf9ecac7fc494cef771d2b6b093da6531  runs/direct-search/frontier-portal-exploitation-20260729/connector32-de764-df0b/aggregate-report.json
```

Pairwise alignment aggregates:

```text
a3a0e1c7799b53b07d3dfea80528529e234bb85dabaf20bae22339cb65413c38  runs/direct-search/frontier-factor-class-expansion-20260728/eb138a-bridges/gram-aut-exhaustive/aggregate-report.json
9e572ab1503d67029170ec24b7327be0430a92237f09f4ce504b7f6b0d5e0a77  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-1e4b/aggregate-report.json
50c971099d7f1d2514b6565614446ed6e3cad776c6fdc419e4f78858283b4619  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-4072/aggregate-report.json
5eedc304d9728a9d06cfae6cec2aed1256ea7357b8dbae646584c193092d3960  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-9035/aggregate-report.json
12cf6e50c437b8078e31272324b1ee4b551f2ef84ec1c4d882b9e455b9a5bf61  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-b584/aggregate-report.json
f6197eb06ee025083d864e33450f3a31e220b8311f2d33b2f16e8e971c1ab9ec  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-b64c/aggregate-report.json
fc11232f3cb3d32704b990e8b11cecc7c7e62506ae75db82a52f28d02357c05d  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-db2c/aggregate-report.json
182765a9e32af28693acdcf870f45864b78f8af56a93230b335012deb736bc13  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-de764/aggregate-report.json
fa797b3849608c2024b6b6f6d4c51cf07f49df68e41402cb69d1e8e0a8a734ba  runs/direct-search/frontier-portal-exploitation-20260729/exact-from-df0b/aggregate-report.json
```

Connector replay/provenance sidecars:

```text
7c6706dc4af26ebd1e820e01bd2fde7eb7f884f294cc52fee6c3fdef9dd9ac30  runs/direct-search/frontier-portal-exploitation-20260729/connector32-de764-1e4b/provenance.json
318df7afe60853a97fa08c503c580b32fff199510fb7988ef143a1492acf84c2  runs/direct-search/frontier-portal-exploitation-20260729/connector32-1e4b-de764-alt/provenance.json
2118359fdee1fb690abfba36dfd9cb55f010ac79c8bc7953fe5703033fe5d1db  runs/direct-search/frontier-portal-exploitation-20260729/connector32-1e4b-9035/provenance.json
c91f546e8cbe49282df7917f7223d0299373fe5908c5b39e4141b19497d8a436  runs/direct-search/frontier-portal-exploitation-20260729/connector32-de764-df0b/provenance.json
```

Completed stochastic logs and retained verifier records:

```text
6bbd22ae7d3c936d36d75fd7c0e229de6b2155aa3a869f579ae1b5a6531644e9  runs/direct-search/frontier-new-portals-core-20260729/de764-arm0/run.jsonl
9af7d4526dcdae60def4f040aa6a923d4d62978ef2b7b518c6295f5833101ef9  runs/direct-search/frontier-new-portals-core-20260729/de764-arm0/best.receipt.json
9efc9191ce7d6631877f66c824a3cba1caf57a72da8fb866e41785350de1d89a  runs/direct-search/frontier-new-portals-core-20260729/1e4b-arm0/run.jsonl
3e4e777581f1f1ff63d66a5d8659097cfa9cbd505ef08f30c0cfd7cb18b333ab  runs/direct-search/frontier-new-portals-core-20260729/1e4b-arm0/best.receipt.json
08125549719a1b106e41038554e5468c461860f74dc6e4943f232535f61b8c4e  runs/direct-search/frontier-portal-fiber-10class-20260729/seed7292303/run.jsonl
6e261200abd783fa0e1bf022c01ff3d389a8724340b25b12db454093f29e7b30  runs/direct-search/frontier-portal-fiber-10class-20260729/seed7292303/checkpoint.json
f2593ddfe3a9c78b9da9f9be79fc9942f5bf62b395a52d089a1faf42638b0fea  runs/direct-search/frontier-portal-fiber-10class-20260729/seed7292303/best.receipt.json
d366aad49e6a224bb18b235fa9fc6b56258027600c4755ead9c4574107c50294  runs/direct-search/frontier-portal-fiber-10class-20260729/seed7292303/first-tie.receipt.json
e73a2f658c2b923e591bfaefbcb9bb396926372f3d145ed724b6d751c61d889c  runs/direct-search/frontier-portal-fiber-10class-20260729/seed7292303/baseline-plus-first-tie-audit.json
```

The distinction between complete exact searches and finite heuristic runs is
essential: the radius-three and four 32-cube statements are exhaustive only
in their explicitly pinned domains; the polar/DR, core, and neutral-fiber
results are observed finite-run outcomes.

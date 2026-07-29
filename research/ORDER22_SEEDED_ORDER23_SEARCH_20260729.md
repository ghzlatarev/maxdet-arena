# Order-22-seeded order-23 search

Three four-arm waves tested exact order-22 border representatives as order-23
search seeds. Each arm ran for 900 seconds. All 12 terminal matrices passed
`arena verify`; none exceeded the effective frontier
`2,779,447,296,000,000`.

## Results

| wave | arms | best `|det|` | current H/HT classes | new versus prior waves | pair assignments | cross assignments |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| component3 | 4 | `2,694,905,856,000,000` | 3 | 3 | 6,366,953,472 | 3,328,180,224 |
| component4 | 4 | `2,694,905,856,000,000` | 4 | 3 | 6,366,953,472 | 3,328,180,224 |
| component5 | 4 | `2,687,827,968,000,000` | 4 | 2 | 4,244,635,648 | 2,218,786,816 |
| **total** | **12** | **`2,694,905,856,000,000`** | — | **8** | **16,978,542,592** | **8,875,147,264** |

The best gap was `84,541,440,000,000`. The waves consumed 10,800 scheduled
arm-seconds. Equal scores did not imply equivalent matrices: component4
produced two distinct H/HT classes at its best score, and component5 produced
three distinct classes at its best score.

Classification used the pinned `pynauty==2.8.8.1` certificate over all
`23^2` dephased pivots, with transpose included in the HT certificate.
Cross-wave comparison removed one already-settled component4 output and two
already-settled component5 outputs. Every one of the eight genuinely new
H/HT representatives received one complete pair pass and one complete
row-column cross pass. The 16 passes covered exactly **25,853,689,856**
assignments and found no strict improvement.

## Frozen manifests

| wave | manifest | SHA-256 |
| --- | --- | --- |
| component3 | `runs/direct-search/order22-component3-order23-search-20260729/manifest.json` | `dd4fbf17ca825fb51a707a637501ced18f00939588abb7c00138f09bee570a98` |
| component4 | `runs/direct-search/order22-component4-order23-search-20260729/manifest.json` | `6e9828894af5167837b5be3a7e4919885369043e002c33d9aadc8d8ef6cad2a5` |
| component5 | `runs/direct-search/order22-component5-order23-search-20260729/manifest.json` | `0925e6a8caa4d40e3b4c6b0f505aecf3d845e2a7e8fe142845fb76cb14d9df48` |

Each manifest pins the starting matrices, search and refinement binaries,
source hashes, terminal JSONL logs, exact receipts, H/HT reports, pass counts,
and cross-wave skip decisions.

## Claim boundary

The 900-second arms are stochastic and do not exhaust their basins. Pair and
cross passes are complete only for their stated local move families around the
retained representatives. H/HT deduplication is exact for the supplied
outputs, not a classification of all order-23 matrices. Arena receipts verify
matrix scores only; these negative results prove neither global optimality nor
a world record.

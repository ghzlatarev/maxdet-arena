# Exact connector audit with repeated block pairs

**Run date:** 2026-07-28

This extends `GRAM_CONNECTOR12_20260728.md` by allowing the three `K2,2`
connectors to reuse a base-block pair. Their edge sets must remain pairwise
disjoint, so every retained graph still has exactly 12 added edges.

## Exact labeled graph count

The edge-count distribution across occupied block pairs is one of:

```text
4 + 4 + 4
8 + 4
12.
```

Deduplicating local unions by their exact bipartite edge mask gives:

| Block-pair type | Two-connector 8-edge unions | Three-connector 12-edge unions |
|---|---:|---:|
| `K3`--`K4` | 21 | 0 |
| `K4`--`K4` | 174 | 36 |

There are 450 labeled single connectors: five `K3`--`K4` pairs with 18
connectors each and ten `K4`--`K4` pairs with 36 each. Therefore:

```text
distinct pairs:
  12,072,240

one repeated pair:
  5 * 21  * (450 - 18)
+ 10 * 174 * (450 - 36)
= 765,720

one tripled pair:
  5 * 0 + 10 * 36
= 360

total unique labeled graphs:
  12,072,240 + 765,720 + 360
= 12,838,320.
```

The three edge-count types are disjoint. In the `8+4` type, the pair carrying
eight edges is uniquely identifiable. Within that pair, exact edge-mask
deduplication removes every alternate `K2,2` decomposition before the graph is
counted.

## Exact symmetry quotient

For each graph, the engine minimizes the established connector canonical key
over every `K2,2` decomposition of its repeated-pair relation. The key already
minimizes over connector order, endpoint orientation, permutations inside each
base block, and permutations of the five `K4` blocks. Taking the least key
over all decompositions is therefore a complete graph-orbit invariant:

- equivalent graphs map their full decomposition sets bijectively and have
  the same minimum;
- equal minima exhibit equivalent connector decompositions, whose edge unions
  give an explicit graph equivalence.

Under

```text
S3 x (S4 wr S5)
```

the 12,838,320 unique labeled graphs reduce to exactly:

```text
113 orbits.
```

This consists of the previous 73 distinct-pair orbits plus 40 new orbits. The
retained orbit multiplicities sum to 12,838,320, and every orbit size divides
the full group order `5,733,089,280`.

## Complete determinant result

Every one of the 113 representatives has a positive exact determinant. The
largest is:

```text
8,296,580,272,619,520,000,000,000,000,000,
```

which is not a square. Complete statistics are:

```text
positive determinants              113
determinants above frontier^2        29
perfect-square determinants           1
frontier square ties                   1
perfect squares above frontier         0
qualified factor-routing survivors     0
```

The sole square is still the published distinct-pair orbit:

```text
7,725,327,271,241,711,616,000,000,000,000
= 2,779,447,296,000,000^2.
```

No repeated-pair orbit reaches the Hasse, sign-shell, or exact-factor stages.
The generated compatibility snapshot and both downstream reports are complete
empty selections rather than timeouts.

## Reproduce

```sh
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wshadow -Werror \
  research/gram_connector12_reuse.cpp \
  -o /tmp/gram_connector12_reuse

/tmp/gram_connector12_reuse \
  --output \
    runs/direct-search/gram-connector12-reuse-complete-20260728.json \
  --route-snapshot \
    runs/direct-search/gram-connector12-reuse-complete-20260728-route.json

python3 research/gram_hasse.py \
  --snapshot \
    runs/direct-search/gram-connector12-reuse-complete-20260728-route.json \
  --all-hits \
  --output \
    runs/direct-search/gram-connector12-reuse-complete-20260728-hasse.json

/tmp/gram_shell_filter \
  --snapshot \
    runs/direct-search/gram-connector12-reuse-complete-20260728-route.json \
  --all-hits \
  --output \
    runs/direct-search/gram-connector12-reuse-complete-20260728-shell.json
```

The primary report SHA-256 is:

```text
08e0260a6169f20e806b083fd12b7ac2a658d665180ccc8ab56e535a68aff57b
```

The optimized complete run took 9.358 seconds. A full ASan/UBSan replay
finished with the same 12,838,320 / 113 / 1 / 0 counts and a semantically
identical report.

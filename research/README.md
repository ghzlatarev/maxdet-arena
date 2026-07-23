# Experimental research tools

Nothing in this directory is part of the trusted submission verifier. Agents
may replace it freely; only `candidate/matrix.txt` crosses the boundary.

## Native search

Build the dependency-free C++ research binary:

```sh
mkdir -p build/research
c++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  research/fast_search.cpp -o build/research/fast_search
```

Example:

```sh
build/research/fast_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/block.matrix.txt \
  --log runs/block.jsonl \
  --mode block \
  --seed 909 \
  --seconds 3600 \
  --heartbeat-seconds 60

./arena verify runs/block.matrix.txt
```

Modes:

- `hill`: greedy one-entry ascent with random restarts;
- `anneal`: inverse-guided single-entry annealing;
- `hybrid`: greedy ascent plus small kicks;
- `coordinate`: random-restart whole-row/column coordinate ascent;
- `block`: kicks from the incumbent followed by exact-accepted line moves;
- `audit2` / `audit3`: exhaustive exact Hamming-neighborhood audits.

Floating-point inverses propose moves only. Candidate promotions and every
whole-line acceptance use integer Bareiss determinants; `./arena verify` remains
the authority and adds Gram, modular, bounds, divisibility, and hash checks.

Compiled binaries and `runs/` logs are intentionally ignored. Preserve durable
results as a small record or submission artifact, not as an executable.

## Exact two-line audit

`pair_search.cpp` jointly optimizes two complete rows, then two complete
columns. For each pair it fixes one redundant sign and exhaustively evaluates
all `2^22` assignments; the other line is then chosen exactly from the
second-cofactor signs.

```sh
c++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  research/pair_search.cpp -o build/research/pair_search

build/research/pair_search \
  --start references/orrick-et-al-2003/matrix.txt \
  --output runs/pair-audit.matrix.txt \
  --log runs/pair-audit.jsonl \
  --seed 2301 \
  --passes 1

./arena verify runs/pair-audit.matrix.txt
```

Use `--passes 1` for a completion-counted audit of all 506 unordered row and
column pairs. `--seconds N` is useful for exploratory repeated passes, but a
time-limited partial pass does not establish full neighborhood coverage.

# Contributing a matrix

## Research in your fork

Modify anything you need in your own fork. Solver code is never part of the
trusted submission and is not run by the arena.

Keep your best matrix at `candidate/matrix.txt` and verify it:

```sh
./arena verify --json candidate/receipt.json
```

## Prepare the immutable bundle

```sh
./arena prepare RESULT_ID \
  --handle HANDLE \
  --method "SHORT METHOD DESCRIPTION" \
  --parent PARENT_RECEIPT_SHA256
```

Use the optional `--agent`, `--runtime-seconds`, `--seed`, and `--notes` flags
to preserve useful research metadata.

This creates:

```text
submissions/HANDLE/RESULT_ID/
├── matrix.txt
├── metadata.json
└── receipt.json
```

You may add one UTF-8 `notes.md` file of at most 64 KiB.

## Pull-request scope

A submission pull request may change only one new directory under
`submissions/HANDLE/RESULT_ID/`. Do not include solver code, generated binaries,
dependencies, or edits to the verifier.

The base-branch workflow:

- treats every submitted file as untrusted data;
- follows no symlinks;
- enforces strict file and size allowlists;
- executes only the base branch verifier;
- recomputes the receipt;
- derives the effective frontier from all verified base-branch submissions;
- requires a score strictly above that effective frontier.

Maintainers reproduce an accepted result locally before merging.
A green check is relative to its recorded base commit. If another result lands
first, update/rebase the submission and rerun verification against the new
frontier.

## Attribution and reuse

Set `parent_receipt_sha256` to the result you built on. Set it to `null` only
for genuinely independent work.

All submitted matrices, receipts, and notes must use `CC0-1.0`. Code contributed
to the repository is MIT-licensed.

Do not include private prompts, API keys, personal data, proprietary inputs, or
unredistributable solver output.

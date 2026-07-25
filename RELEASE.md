# Private-to-public release gate

Do not call the arena open until every item below is complete.

## Scientific and harness gate

- [x] Multi-strategy private dogfood has run for at least two hours.
- [x] Every retained matrix passes the Python exact verifier.
- [x] `npm ci && npm run check` passes from a clean clone.
- [x] The data-only pull-request integration tests accept an improvement and
      reject a tie, code edit, malformed bundle, and symlink.
- [x] The published comparison point has a final literature/provenance review.
- [x] Dogfood results and known harness limitations are summarized in
      `research/DOGFOOD.md`.

## Repository gate

- [x] Create the intended GitHub repository and push the exact reviewed commit.
- [x] Make the repository public after explicit release approval.
- [x] Keep Actions permissions read-only by default.
- [ ] Protect `main`; block force pushes and require the exact-verification
      checks for public submissions.
- [ ] Require submission branches to be current with `main` before merge so
      two simultaneously green results cannot race an older frontier.
- [x] Enable Discussions or remove the Discussions link.
- [x] Verify the one-line clone/Codex command against the real URL.

## Site gate

- [ ] Change `data/frontier.json` status from `private-dogfooding` to `open`.
- [x] Push the same reviewed commit to the connected Sites source repository.
- [x] Save and inspect a private version before changing access.
- [x] Make any public deployment only with explicit release approval.

## Current public-preview state

- The GitHub repository is public and keeps read-only Actions permissions by
  default. Protect `main` before calling the arena open for submissions.
- GitHub Pages deploys the static site only after the exact repository tests,
  type-check, and production export pass. The connected Sites project remains
  owner-only and is not deployed.
- `data/frontier.json` remains `private-dogfooding`. This is a public preview,
  not yet an open-arena claim.

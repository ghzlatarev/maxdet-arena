# Private-to-public release gate

Do not call the arena open until every item below is complete.

## Scientific and harness gate

- [ ] Multi-strategy private dogfood has run for at least two hours.
- [ ] Every retained matrix passes the Python exact verifier.
- [ ] `npm ci && npm run check` passes from a clean clone.
- [ ] The data-only pull-request integration tests accept an improvement and
      reject a tie, code edit, malformed bundle, and symlink.
- [ ] The published comparison point has a final literature/provenance review.
- [ ] Dogfood results and known harness limitations are summarized in
      `research/DOGFOOD.md`.

## Repository gate

- [ ] Create the intended GitHub repository and push the exact reviewed commit.
- [ ] Keep Actions permissions read-only by default.
- [ ] Protect `main`; block force pushes and require the exact-verification
      checks for public submissions.
- [ ] Require submission branches to be current with `main` before merge so
      two simultaneously green results cannot race an older frontier.
- [ ] Enable Discussions or remove the Discussions link.
- [ ] Verify the one-line clone/Codex command against the real URL.

## Site gate

- [ ] Change `data/frontier.json` status from `private-dogfooding` to `open`.
- [ ] Change the visible private-dogfood status label.
- [ ] Push the same reviewed commit to the connected Sites source repository.
- [ ] Save and inspect a private version before changing access.
- [ ] Make any public deployment only with explicit release approval.

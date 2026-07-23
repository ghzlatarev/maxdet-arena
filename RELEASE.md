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

- [ ] Create the intended GitHub repository and push the exact reviewed commit.
- [x] Keep Actions permissions read-only by default.
- [ ] Protect `main`; block force pushes and require the exact-verification
      checks for public submissions.
- [ ] Require submission branches to be current with `main` before merge so
      two simultaneously green results cannot race an older frontier.
- [x] Enable Discussions or remove the Discussions link.
- [ ] Verify the one-line clone/Codex command against the real URL.

## Site gate

- [ ] Change `data/frontier.json` status from `private-dogfooding` to `open`.
- [x] Push the same reviewed commit to the connected Sites source repository.
- [x] Save and inspect a private version before changing access.
- [ ] Make any public deployment only with explicit release approval.

## Current private state

- The intended GitHub repository exists with private visibility and read-only
  Actions defaults. Its first source push is pending a GitHub OAuth token with
  `workflow` scope; do not omit the pinned verification workflow to bypass it.
- The connected Sites project is owner-only. A saved source version is not a
  deployment and does not create a public URL.
- `data/frontier.json` remains `private-dogfooding`. Public access, an open
  frontier status, and deployment all require explicit release approval.

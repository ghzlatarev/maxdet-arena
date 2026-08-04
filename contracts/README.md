# MaxDet order-23 Sepolia bounty

This is a one-shot, ownerless testnet pilot. Anyone may irrevocably donate
Sepolia ETH. The first successful prior-block committer to reveal a 23 × 23
sign matrix with exact absolute determinant at least
`2,779,447,300,194,304` earns the entire contract balance.

The contract has no multisig, oracle, admin, upgrade, timeout, refund, token, or
yield strategy. Testnet ETH has no financial value.

## One acceptance rule

`claim` accepts 23 unsigned integers. Integer `i` encodes matrix row `i`; bit
`j` is `+1` when set and `-1` when clear. Bits above 22 are rejected. The
contract sign-normalizes the matrix, reduces it to a 22 × 22 binary core, and
uses fraction-free Bareiss elimination to compute

```text
abs(det(A)) = 2^22 * abs(det(binary core))
```

The published comparison matrix scores `2,779,447,296,000,000`, exactly one
determinant lattice step below the minimum payout score, so a tie cannot claim.

## Test

Install [Foundry](https://getfoundry.sh/) and run:

```sh
forge fmt --check && forge test --gas-report && forge build --sizes
```

The production-compiler gate is 10 million gas for a claim. The measured
reference-vector claim path is under 4 million gas, and runtime bytecode is
about 3.7 kB.

## Deploy to Sepolia

Use a dedicated encrypted testnet keystore, never a plaintext key:

```sh
cast wallet import mathfast-sepolia --interactive
```

Then deploy the immutable source:

```sh
forge create contracts/src/MaxDetBounty23.sol:MaxDetBounty23 --chain 11155111 --rpc-url "$SEPOLIA_RPC_URL" --account mathfast-sepolia --broadcast --verify --verifier sourcify
```

Before accepting donations, update `contracts/deployments/sepolia.json` with
the address, transaction, deployment block, and `live` status. Verify the source
and compiler settings on an explorer, compare `keccak256` of the deployed
runtime bytecode with `expected_runtime_codehash`, run the registry tests, and
confirm the site enables donations only after its live code check passes.

## Claim before publishing

First verify and encode the private winning matrix locally:

```sh
./arena verify candidate/matrix.txt
MATRIX=$(python3 tools/encode_onchain_matrix.py candidate/matrix.txt --require-winning --cast-array)
```

Keep the matrix and salt private until the commitment is mined. Compute the
commitment locally; do not call a third-party RPC helper with an unrevealed
matrix:

```sh
CLAIMANT=$(cast wallet address --account "$MAXDET_SEPOLIA_ACCOUNT")
RECIPIENT=${MAXDET_PAYOUT_RECIPIENT:-$CLAIMANT}
SALT=0x$(openssl rand -hex 32)
ENCODED=$(cast abi-encode "commitment(uint256,address,address,address,uint32[23],bytes32)" 11155111 "$MAXDET_BOUNTY_ADDRESS" "$CLAIMANT" "$RECIPIENT" "$MATRIX" "$SALT")
COMMITMENT=$(cast keccak "$ENCODED")
cast send "$MAXDET_BOUNTY_ADDRESS" "commitClaim(bytes32)" "$COMMITMENT" --chain 11155111 --rpc-url "$SEPOLIA_RPC_URL" --account "$MAXDET_SEPOLIA_ACCOUNT"
```

Wait for the commitment receipt and preferably at least two confirmations.
Then reveal through trusted RPC endpoints. The contract recomputes the score
and attempts to pay the committed recipient immediately:

```sh
cast send "$MAXDET_BOUNTY_ADDRESS" "claim(uint32[23],bytes32,address)" "$MATRIX" "$SALT" "$RECIPIENT" --chain 11155111 --rpc-url "$SEPOLIA_RPC_URL" --account "$MAXDET_SEPOLIA_ACCOUNT"
```

If that recipient rejects ETH, the winning result remains final and the prize
becomes a winner-controlled credit. Redirect it with:

```sh
RECOVERY_RECIPIENT=0xYOUR_SAFE_PAYOUT_ADDRESS
cast send "$MAXDET_BOUNTY_ADDRESS" "withdrawPrize(address)" "$RECOVERY_RECIPIENT" --chain 11155111 --rpc-url "$SEPOLIA_RPC_URL" --account "$MAXDET_SEPOLIA_ACCOUNT"
```

Only after the successful reveal should the matrix be pushed to GitHub.

## Donation semantics

- Donations are gifts and can never be withdrawn or refunded.
- The bounty must have a nonzero balance before it can be solved.
- The claimant and payout recipient may differ, but both are fixed by the
  commitment.
- A rejected automatic payout cannot reopen the problem; only the recorded
  winner can redirect the deferred prize.
- Normal donations after the first valid claim revert.
- If nobody meets the fixed threshold, the funds remain locked indefinitely.
- ETH forcibly sent after closure can remain stranded.

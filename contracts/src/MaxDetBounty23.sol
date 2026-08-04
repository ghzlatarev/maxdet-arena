// SPDX-License-Identifier: MIT
pragma solidity 0.8.36;

/// @title MaxDetBounty23
/// @notice A one-shot, donation-funded bounty for an order-23 MaxDet improvement.
/// @dev Each uint32 encodes one matrix row. Bit j is +1 when set and -1 when
///      clear. Only the low 23 bits may be used. There is no owner, timeout,
///      refund, upgrade, oracle, or discretionary verifier.
contract MaxDetBounty23 {
    struct Commitment {
        bytes32 digest;
        uint64 blockNumber;
    }

    uint256 public constant ORDER = 23;
    uint256 public constant ROW_MASK = (2 ** ORDER) - 1;
    uint256 public constant DETERMINANT_LATTICE_STEP = 2 ** (ORDER - 1);
    uint256 public constant PUBLISHED_FRONTIER = 2_779_447_296_000_000;
    uint256 public constant MINIMUM_WINNING_DETERMINANT = PUBLISHED_FRONTIER + DETERMINANT_LATTICE_STEP;
    uint256 public constant AUTOMATIC_PAYOUT_GAS = 30_000;
    string public constant CHALLENGE_ID = "maxdet-23-v1";

    mapping(address claimant => Commitment commitment) public commitments;

    bool public solved;
    address public winner;
    address public payoutRecipient;
    uint256 public winningDeterminant;
    bytes32 public winningMatrixHash;
    uint256 public totalDonated;
    uint256 public claimablePrize;

    error AlreadySolved();
    error BelowMinimum(uint256 determinant, uint256 minimum);
    error CommitmentMismatch(bytes32 committed, bytes32 expected);
    error CommitmentTooFresh(uint256 committedAt, uint256 currentBlock);
    error EmptyDonation();
    error EmptyPrize();
    error InvalidRow(uint256 rowIndex, uint256 encodedRow);
    error InvalidRecipient();
    error NotWinner();
    error NonExactDivision();
    error PayoutFailed();

    event ClaimCommitted(address indexed claimant, bytes32 indexed commitment);
    event Donated(address indexed donor, uint256 amount, uint256 prizeBalance);
    event PrizeClaimed(
        address indexed claimant,
        address indexed recipient,
        bytes32 indexed matrixHash,
        uint256 determinant,
        uint256 payout
    );
    event PrizeWithdrawn(address indexed winner, address indexed recipient, uint256 payout);

    receive() external payable {
        donate();
    }

    /// @notice Irrevocably add ETH to the one-shot prize.
    function donate() public payable {
        if (solved) revert AlreadySolved();
        if (msg.value == 0) revert EmptyDonation();

        totalDonated += msg.value;
        emit Donated(msg.sender, msg.value, address(this).balance);
    }

    /// @notice Record a hidden claim at least one block before revealing it.
    /// @dev The delay prevents a copied reveal transaction from stealing payout.
    function commitClaim(bytes32 commitment) external {
        if (solved) revert AlreadySolved();

        commitments[msg.sender] = Commitment(commitment, uint64(block.number));
        emit ClaimCommitted(msg.sender, commitment);
    }

    /// @notice Verify a winner and attempt to pay the full prize to its committed recipient.
    function claim(uint32[23] calldata rows, bytes32 salt, address payable recipient) external {
        if (solved) revert AlreadySolved();
        if (recipient == address(0)) revert InvalidRecipient();

        bytes32 commitment = _commitmentFor(msg.sender, recipient, rows, salt);
        Commitment memory pending = commitments[msg.sender];
        if (pending.digest != commitment) {
            revert CommitmentMismatch(pending.digest, commitment);
        }
        uint256 committedAt = pending.blockNumber;
        if (block.number <= committedAt) {
            revert CommitmentTooFresh(committedAt, block.number);
        }

        uint256 determinant = scoreMatrix(rows);
        uint256 minimum = minimumWinningDeterminant();
        if (determinant < minimum) {
            revert BelowMinimum(determinant, minimum);
        }

        bytes32 matrixHash = keccak256(abi.encode(rows));
        uint256 payout = address(this).balance;
        if (payout == 0) revert EmptyPrize();

        delete commitments[msg.sender];
        solved = true;
        winner = msg.sender;
        payoutRecipient = recipient;
        winningDeterminant = determinant;
        winningMatrixHash = matrixHash;

        // Record the full credit before interacting with recipient code. A
        // bounded best-effort push keeps ordinary EOA claims atomic; failed
        // delivery leaves the winner-controlled credit intact.
        claimablePrize = payout;
        emit PrizeClaimed(msg.sender, recipient, matrixHash, determinant, payout);
        if (_tryAutomaticPayout(recipient, payout)) claimablePrize = 0;
    }

    /// @notice Redirect a deferred prize after the committed recipient rejected ETH.
    function withdrawPrize(address payable recipient) external {
        if (msg.sender != winner) revert NotWinner();
        if (recipient == address(0)) revert InvalidRecipient();

        uint256 payout = claimablePrize;
        if (payout == 0) revert EmptyPrize();

        claimablePrize = 0;
        payoutRecipient = recipient;
        emit PrizeWithdrawn(msg.sender, recipient, payout);

        (bool paid,) = recipient.call{value: payout}("");
        if (!paid) revert PayoutFailed();
    }

    /// @dev Uses no return-data buffer, so a recipient cannot return-bomb the
    ///      verifier. The fixed budget is enough for an EOA and deliberately
    ///      treats more complex recipients as deferred withdrawals.
    function _tryAutomaticPayout(address recipient, uint256 payout) private returns (bool paid) {
        uint256 payoutGas = AUTOMATIC_PAYOUT_GAS;
        assembly ("memory-safe") {
            paid := call(payoutGas, recipient, payout, 0, 0, 0, 0)
        }
    }

    /// @notice Return the exact absolute determinant represented by `rows`.
    function scoreMatrix(uint32[23] calldata rows) public pure returns (uint256) {
        int256 reduced = _reducedDeterminant(rows);
        uint256 absoluteReduced = uint256(reduced < 0 ? -reduced : reduced);
        return absoluteReduced * DETERMINANT_LATTICE_STEP;
    }

    /// @notice Reproduce the commitment hash locally before calling commitClaim.
    function commitmentFor(address claimant, address recipient, uint32[23] calldata rows, bytes32 salt)
        external
        view
        returns (bytes32)
    {
        return _commitmentFor(claimant, recipient, rows, salt);
    }

    /// @notice The only mathematical acceptance threshold.
    function minimumWinningDeterminant() public pure virtual returns (uint256) {
        return MINIMUM_WINNING_DETERMINANT;
    }

    function _commitmentFor(address claimant, address recipient, uint32[23] calldata rows, bytes32 salt)
        internal
        view
        returns (bytes32)
    {
        return keccak256(abi.encode(block.chainid, address(this), claimant, recipient, rows, salt));
    }

    /// @dev Sign-normalize the first row and column, subtract the first row,
    ///      and factor 2 from each remaining row. This reduces the problem to
    ///      an exact 22 x 22 binary determinant:
    ///
    ///          abs(det(A)) = 2^22 * abs(det(B)).
    function _reducedDeterminant(uint32[23] calldata rows) internal pure returns (int256) {
        for (uint256 i; i < ORDER;) {
            if (uint256(rows[i]) > ROW_MASK) {
                revert InvalidRow(i, rows[i]);
            }
            unchecked {
                ++i;
            }
        }

        // Solidity memory is zero-initialized; parity-one cells are filled below.
        int256[484] memory matrix;
        uint32 firstRow = rows[0];
        uint32 firstCorner = firstRow & 1;

        for (uint256 i = 1; i < ORDER;) {
            uint32 row = rows[i];
            uint32 firstColumn = row & 1;
            for (uint256 j = 1; j < ORDER;) {
                uint32 parity = ((row >> j) ^ (firstRow >> j) ^ firstColumn ^ firstCorner) & 1;
                matrix[(i - 1) * 22 + (j - 1)] = int256(uint256(parity));
                unchecked {
                    ++j;
                }
            }
            unchecked {
                ++i;
            }
        }

        return _bareiss22(matrix);
    }

    /// @dev Fraction-free Bareiss elimination with deterministic row pivoting.
    ///      Hadamard bounds put every stored 0/1 minor below 2^31 and every
    ///      multiply-subtract below 2^63 here. Checked int256 arithmetic still
    ///      provides a hard failure if that invariant is ever violated.
    function _bareiss22(int256[484] memory matrix) internal pure returns (int256) {
        int256 previousPivot = 1;
        int256 sign = 1;

        for (uint256 k; k < 21;) {
            uint256 pivotRow = k;
            while (pivotRow < 22 && matrix[pivotRow * 22 + k] == 0) {
                unchecked {
                    ++pivotRow;
                }
            }
            if (pivotRow == 22) return 0;

            if (pivotRow != k) {
                for (uint256 j = k; j < 22;) {
                    (matrix[k * 22 + j], matrix[pivotRow * 22 + j]) = (matrix[pivotRow * 22 + j], matrix[k * 22 + j]);
                    unchecked {
                        ++j;
                    }
                }
                sign = -sign;
            }

            int256 pivot = matrix[k * 22 + k];
            for (uint256 i = k + 1; i < 22;) {
                int256 lower = matrix[i * 22 + k];
                for (uint256 j = k + 1; j < 22;) {
                    int256 numerator = matrix[i * 22 + j] * pivot - lower * matrix[k * 22 + j];
                    if (k != 0) {
                        if (numerator % previousPivot != 0) {
                            revert NonExactDivision();
                        }
                        numerator /= previousPivot;
                    }
                    matrix[i * 22 + j] = numerator;
                    unchecked {
                        ++j;
                    }
                }
                matrix[i * 22 + k] = 0;
                unchecked {
                    ++i;
                }
            }

            previousPivot = pivot;
            unchecked {
                ++k;
            }
        }

        return sign * matrix[483];
    }
}

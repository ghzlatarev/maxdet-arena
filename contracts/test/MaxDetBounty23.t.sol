// SPDX-License-Identifier: MIT
pragma solidity 0.8.36;

import {MaxDetBounty23} from "../src/MaxDetBounty23.sol";

interface Vm {
    function deal(address account, uint256 balance) external;
    function expectRevert(bytes calldata revertData) external;
    function expectRevert(bytes4 selector) external;
    function prank(address caller) external;
    function roll(uint256 newHeight) external;
}

contract MaxDetBounty23Harness is MaxDetBounty23 {
    function minimumWinningDeterminant() public pure override returns (uint256) {
        return PUBLISHED_FRONTIER;
    }
}

contract ReentrantClaimant {
    MaxDetBounty23 private immutable bounty;
    bool public reentrySucceeded;

    constructor(MaxDetBounty23 bounty_) {
        bounty = bounty_;
    }

    function commit(uint32[23] calldata rows, bytes32 salt) external {
        bounty.commitClaim(bounty.commitmentFor(address(this), address(this), rows, salt));
    }

    function reveal(uint32[23] calldata rows, bytes32 salt) external {
        bounty.claim(rows, salt, payable(address(this)));
    }

    receive() external payable {
        (reentrySucceeded,) = address(bounty).call(abi.encodeCall(MaxDetBounty23.commitClaim, (bytes32(uint256(1)))));
    }
}

contract RevertingRecipient {
    receive() external payable {
        revert("reject ETH");
    }
}

contract MaxDetBounty23Test {
    Vm private constant vm = Vm(address(uint160(uint256(keccak256("hevm cheat code")))));

    MaxDetBounty23 private bounty;
    address private constant DONOR = address(0xD0A0);
    address private constant CLAIMANT = address(0xC1A1);
    address payable private constant RECIPIENT = payable(address(0xB0B));
    bytes32 private constant SALT = keccak256("maxdet-test-salt");

    function setUp() public {
        bounty = new MaxDetBounty23();
    }

    function testReferenceMatrixScoresExactly() public view {
        uint32[23] memory rows = _referenceRows();
        _assertEq(bounty.scoreMatrix(rows), 2_779_447_296_000_000);
    }

    function testReferenceMatrixIsOneLatticeStepBelowWinningMinimum() public view {
        _assertEq(bounty.minimumWinningDeterminant() - bounty.scoreMatrix(_referenceRows()), 1 << 22);
    }

    function testAllPositiveMatrixIsSingular() public view {
        uint32[23] memory rows;
        for (uint256 i; i < 23; ++i) {
            rows[i] = uint32((1 << 23) - 1);
        }
        _assertEq(bounty.scoreMatrix(rows), 0);
    }

    function testTwoIdentityMinusOnesMatrixScoresExactly() public view {
        uint32[23] memory rows;
        for (uint256 i; i < 23; ++i) {
            // casting to 'uint32' is safe because i is less than 23.
            // forge-lint: disable-next-line(unsafe-typecast)
            rows[i] = uint32(2 ** i);
        }
        _assertEq(bounty.scoreMatrix(rows), 21 * (1 << 22));
    }

    function testReversePermutationCoreForcesRepeatedPivotSwaps() public view {
        uint32[23] memory rows;
        uint32 mask = uint32((1 << 23) - 1);
        rows[0] = mask;
        for (uint256 coreRow; coreRow < 22; ++coreRow) {
            uint256 reversedColumn = 21 - coreRow;
            rows[coreRow + 1] = mask ^ _bit(reversedColumn + 1);
        }
        _assertEq(bounty.scoreMatrix(rows), 1 << 22);
    }

    function testLateCoreDependencyIsSingular() public view {
        uint32[23] memory rows;
        uint32 mask = uint32((1 << 23) - 1);
        rows[0] = mask;
        for (uint256 coreRow; coreRow < 22; ++coreRow) {
            uint256 coreColumn = coreRow == 21 ? 20 : coreRow;
            rows[coreRow + 1] = mask ^ _bit(coreColumn + 1);
        }
        _assertEq(bounty.scoreMatrix(rows), 0);
    }

    function testPythonDifferentialVectorSeed1() public view {
        _assertEq(bounty.scoreMatrix(_lcgRows(1)), 12_771_655_680);
    }

    function testPythonDifferentialVectorSeed23() public view {
        _assertEq(bounty.scoreMatrix(_lcgRows(23)), 9_302_966_272);
    }

    function testPythonDifferentialVectorSeed2026() public view {
        _assertEq(bounty.scoreMatrix(_lcgRows(2026)), 3_305_111_552);
    }

    function testScoreIsInvariantUnderSignFlipsAndPermutations() public view {
        uint32[23] memory rows = _referenceRows();
        uint256 expected = bounty.scoreMatrix(rows);
        uint32 mask = uint32((1 << 23) - 1);

        rows[4] ^= mask;
        for (uint256 i; i < 23; ++i) {
            rows[i] ^= uint32(1 << 9);
        }
        (rows[2], rows[17]) = (rows[17], rows[2]);
        for (uint256 i; i < 23; ++i) {
            uint32 left = (rows[i] >> 5) & 1;
            uint32 right = (rows[i] >> 19) & 1;
            if (left != right) rows[i] ^= uint32((1 << 5) | (1 << 19));
        }

        _assertEq(bounty.scoreMatrix(rows), expected);
    }

    function testScoreIsInvariantUnderTranspose() public view {
        uint32[23] memory rows = _referenceRows();
        uint32[23] memory transposed;
        for (uint256 row; row < 23; ++row) {
            for (uint256 column; column < 23; ++column) {
                if ((rows[row] & _bit(column)) != 0) {
                    transposed[column] |= _bit(row);
                }
            }
        }
        _assertEq(bounty.scoreMatrix(transposed), 2_779_447_296_000_000);
    }

    function testRejectsBitsOutsideOrder23() public {
        uint32[23] memory rows = _referenceRows();
        rows[7] |= uint32(1 << 23);
        vm.expectRevert(abi.encodeWithSelector(MaxDetBounty23.InvalidRow.selector, 7, rows[7]));
        bounty.scoreMatrix(rows);
    }

    function testAnyoneCanDonate() public {
        vm.deal(DONOR, 3 ether);
        vm.prank(DONOR);
        bounty.donate{value: 2 ether}();

        _assertEq(address(bounty).balance, 2 ether);
        _assertEq(bounty.totalDonated(), 2 ether);
    }

    function testEmptyDonationReverts() public {
        vm.expectRevert(MaxDetBounty23.EmptyDonation.selector);
        bounty.donate();
    }

    function testClaimNeedsACommitment() public {
        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = bounty.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.expectRevert(abi.encodeWithSelector(MaxDetBounty23.CommitmentMismatch.selector, bytes32(0), commitment));
        vm.prank(CLAIMANT);
        bounty.claim(rows, SALT, RECIPIENT);
    }

    function testClaimNeedsANewerBlock() public {
        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = bounty.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        bounty.commitClaim(commitment);

        vm.expectRevert(abi.encodeWithSelector(MaxDetBounty23.CommitmentTooFresh.selector, block.number, block.number));
        vm.prank(CLAIMANT);
        bounty.claim(rows, SALT, RECIPIENT);
    }

    function testClaimRejectsZeroRecipient() public {
        vm.expectRevert(MaxDetBounty23.InvalidRecipient.selector);
        vm.prank(CLAIMANT);
        bounty.claim(_referenceRows(), SALT, payable(address(0)));
    }

    function testPayoutRecipientCannotChangeAfterCommit() public {
        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = bounty.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        bounty.commitClaim(commitment);
        vm.roll(block.number + 1);

        address payable replacement = payable(address(0xCAFE));
        bytes32 replacementCommitment = bounty.commitmentFor(CLAIMANT, replacement, rows, SALT);
        vm.expectRevert(
            abi.encodeWithSelector(MaxDetBounty23.CommitmentMismatch.selector, commitment, replacementCommitment)
        );
        vm.prank(CLAIMANT);
        bounty.claim(rows, SALT, replacement);
    }

    function testCopiedRevealCannotUseAnotherClaimantsCommitment() public {
        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = bounty.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        bounty.commitClaim(commitment);

        address copier = address(0xC0FFEE);
        vm.prank(copier);
        bounty.commitClaim(commitment);
        vm.roll(block.number + 1);

        bytes32 copierCommitment = bounty.commitmentFor(copier, RECIPIENT, rows, SALT);
        vm.expectRevert(
            abi.encodeWithSelector(MaxDetBounty23.CommitmentMismatch.selector, commitment, copierCommitment)
        );
        vm.prank(copier);
        bounty.claim(rows, SALT, RECIPIENT);
    }

    function testPublishedFrontierCannotClaimProductionBounty() public {
        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = bounty.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        bounty.commitClaim(commitment);
        vm.roll(block.number + 1);

        vm.expectRevert(
            abi.encodeWithSelector(MaxDetBounty23.BelowMinimum.selector, 2_779_447_296_000_000, 2_779_447_300_194_304)
        );
        vm.prank(CLAIMANT);
        bounty.claim(rows, SALT, RECIPIENT);
    }

    function testQualifyingMatrixCannotCloseAnEmptyPrize() public {
        MaxDetBounty23Harness harness = new MaxDetBounty23Harness();
        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = harness.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        harness.commitClaim(commitment);
        vm.roll(block.number + 1);

        vm.expectRevert(MaxDetBounty23.EmptyPrize.selector);
        vm.prank(CLAIMANT);
        harness.claim(rows, SALT, RECIPIENT);
    }

    function testValidClaimPaysTheEntirePrizeAndClosesBounty() public {
        MaxDetBounty23Harness harness = new MaxDetBounty23Harness();
        vm.deal(DONOR, 5 ether);
        vm.prank(DONOR);
        harness.donate{value: 5 ether}();

        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = harness.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        harness.commitClaim(commitment);
        vm.roll(block.number + 1);

        uint256 recipientBefore = RECIPIENT.balance;
        vm.prank(CLAIMANT);
        harness.claim(rows, SALT, RECIPIENT);

        _assertTrue(harness.solved());
        _assertEq(harness.winner(), CLAIMANT);
        _assertEq(harness.payoutRecipient(), RECIPIENT);
        _assertEq(harness.winningDeterminant(), 2_779_447_296_000_000);
        _assertEq(harness.claimablePrize(), 0);
        _assertEq(address(harness).balance, 0);
        _assertEq(RECIPIENT.balance, recipientBefore + 5 ether);

        vm.expectRevert(MaxDetBounty23.AlreadySolved.selector);
        harness.donate{value: 1}();

        vm.expectRevert(MaxDetBounty23.AlreadySolved.selector);
        vm.prank(CLAIMANT);
        harness.commitClaim(commitment);

        vm.expectRevert(MaxDetBounty23.AlreadySolved.selector);
        vm.prank(CLAIMANT);
        harness.claim(rows, SALT, RECIPIENT);
    }

    function testClaimExecutionGasStaysBelowTenMillion() public {
        MaxDetBounty23Harness harness = new MaxDetBounty23Harness();
        vm.deal(DONOR, 1 ether);
        vm.prank(DONOR);
        harness.donate{value: 1 ether}();

        uint32[23] memory rows = _referenceRows();
        bytes32 commitment = harness.commitmentFor(CLAIMANT, RECIPIENT, rows, SALT);
        vm.prank(CLAIMANT);
        harness.commitClaim(commitment);
        vm.roll(block.number + 1);

        uint256 gasBefore = gasleft();
        vm.prank(CLAIMANT);
        harness.claim(rows, SALT, RECIPIENT);
        uint256 gasUsed = gasBefore - gasleft();

        require(gasUsed < 10_000_000, "claim gas ceiling exceeded");
    }

    function testPayoutCallbackCannotReenter() public {
        MaxDetBounty23Harness harness = new MaxDetBounty23Harness();
        ReentrantClaimant claimant = new ReentrantClaimant(harness);
        vm.deal(DONOR, 1 ether);
        vm.prank(DONOR);
        harness.donate{value: 1 ether}();

        uint32[23] memory rows = _referenceRows();
        claimant.commit(rows, SALT);
        vm.roll(block.number + 1);
        claimant.reveal(rows, SALT);

        _assertTrue(harness.solved());
        _assertEq(harness.winner(), address(claimant));
        _assertEq(harness.payoutRecipient(), address(claimant));
        _assertEq(address(claimant).balance, 1 ether);
        _assertEq(harness.claimablePrize(), 0);
        require(!claimant.reentrySucceeded(), "reentry unexpectedly succeeded");
    }

    function testRevertingPayoutDefersPrizeAndWinnerCanRedirect() public {
        MaxDetBounty23Harness harness = new MaxDetBounty23Harness();
        RevertingRecipient rejecting = new RevertingRecipient();
        vm.deal(DONOR, 2 ether);
        vm.prank(DONOR);
        harness.donate{value: 2 ether}();

        uint32[23] memory rows = _referenceRows();
        bytes32 rejectedCommitment = harness.commitmentFor(CLAIMANT, address(rejecting), rows, SALT);
        vm.roll(100);
        vm.prank(CLAIMANT);
        harness.commitClaim(rejectedCommitment);
        vm.roll(101);

        vm.prank(CLAIMANT);
        harness.claim(rows, SALT, payable(address(rejecting)));

        _assertTrue(harness.solved());
        _assertEq(harness.winner(), CLAIMANT);
        _assertEq(harness.payoutRecipient(), address(rejecting));
        _assertEq(harness.winningDeterminant(), 2_779_447_296_000_000);
        _assertEq(address(harness).balance, 2 ether);
        _assertEq(harness.claimablePrize(), 2 ether);
        (bytes32 deletedDigest, uint64 deletedBlock) = harness.commitments(CLAIMANT);
        _assertEq(uint256(deletedDigest), 0);
        _assertEq(uint256(deletedBlock), 0);

        vm.expectRevert(MaxDetBounty23.AlreadySolved.selector);
        vm.prank(address(0xC0FFEE));
        harness.commitClaim(rejectedCommitment);

        vm.expectRevert(MaxDetBounty23.NotWinner.selector);
        vm.prank(address(0xC0FFEE));
        harness.withdrawPrize(RECIPIENT);

        vm.expectRevert(MaxDetBounty23.PayoutFailed.selector);
        vm.prank(CLAIMANT);
        harness.withdrawPrize(payable(address(rejecting)));
        _assertEq(harness.claimablePrize(), 2 ether);
        _assertEq(harness.payoutRecipient(), address(rejecting));

        uint256 recipientBefore = RECIPIENT.balance;
        vm.prank(CLAIMANT);
        harness.withdrawPrize(RECIPIENT);
        _assertEq(RECIPIENT.balance, recipientBefore + 2 ether);
        _assertEq(harness.claimablePrize(), 0);
        _assertEq(harness.winner(), CLAIMANT);
        _assertEq(harness.payoutRecipient(), RECIPIENT);
    }

    function _referenceRows() internal pure returns (uint32[23] memory rows) {
        rows = [
            uint32(8_380_805),
            8_362_086,
            8_357_403,
            1_679_244,
            6_692_756,
            3_343_234,
            5_014_914,
            3_443_193,
            4_919_801,
            2_949_180,
            5_406_812,
            2_816_154,
            5_601_562,
            1_980_417,
            6_412_289,
            478_415,
            447_279,
            240_951,
            373_463,
            6_826_831,
            7_353_519,
            5_786_295,
            3_691_863
        ];
    }

    /// @dev Vectors and expected scores are independently generated by the
    /// repository's arbitrary-precision Python Bareiss implementation.
    function _lcgRows(uint256 state) internal pure returns (uint32[23] memory rows) {
        for (uint256 i; i < 23; ++i) {
            state = (1_103_515_245 * state + 12_345) & 0x7fffffff;
            // casting to 'uint32' is safe after the 31-bit mask.
            // forge-lint: disable-next-line(unsafe-typecast)
            rows[i] = uint32(state) & uint32((1 << 23) - 1);
        }
    }

    function _bit(uint256 index) private pure returns (uint32) {
        require(index < 23, "bit index out of range");
        // casting to 'uint32' is safe because index is less than 23.
        // forge-lint: disable-next-line(unsafe-typecast)
        return uint32(2 ** index);
    }

    function _assertEq(uint256 left, uint256 right) private pure {
        require(left == right, "uint mismatch");
    }

    function _assertEq(address left, address right) private pure {
        require(left == right, "address mismatch");
    }

    function _assertTrue(bool value) private pure {
        require(value, "expected true");
    }
}

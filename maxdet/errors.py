"""Typed errors returned by the trusted verifier."""


class ArenaError(Exception):
    """Base class for deterministic user-facing verification errors."""


class ContractError(ArenaError):
    """The challenge contract is missing, malformed, or unsupported."""


class MatrixFormatError(ArenaError):
    """A candidate matrix does not satisfy the strict text format."""


class VerificationError(ArenaError):
    """Independent exact checks disagree or an invariant is violated."""


class SubmissionError(ArenaError):
    """A submission bundle is malformed or unsafe to inspect."""

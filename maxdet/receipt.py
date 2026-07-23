"""Deterministic verification receipts for MaxDet submissions."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .contract import Contract, load_matrix, parse_matrix_bytes
from .errors import VerificationError
from .exact import (
    bareiss_determinant,
    gram_matrix,
    matrix_text,
    modular_determinant,
    normalize_signs,
)

VERIFIER_VERSION = "maxdet-verifier-0.4.0"


@dataclass(frozen=True)
class VerifiedMatrix:
    matrix: list[list[int]]
    receipt: dict[str, Any]

    @property
    def abs_determinant(self) -> int:
        return int(self.receipt["score"]["absolute_determinant"])


def canonical_json_bytes(data: dict[str, Any]) -> bytes:
    return (
        json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        + "\n"
    ).encode("ascii")


def verify_matrix(path: Path, contract: Contract) -> VerifiedMatrix:
    raw_bytes, matrix = load_matrix(path, contract)
    return _verify_parsed_matrix(raw_bytes, matrix, contract)


def verify_matrix_bytes(
    raw_bytes: bytes,
    contract: Contract,
) -> VerifiedMatrix:
    """Verify bytes already obtained through a trusted, no-following read."""

    matrix = parse_matrix_bytes(raw_bytes, contract)
    return _verify_parsed_matrix(raw_bytes, matrix, contract)


def _verify_parsed_matrix(
    raw_bytes: bytes,
    matrix: list[list[int]],
    contract: Contract,
) -> VerifiedMatrix:
    determinant = bareiss_determinant(matrix)
    absolute = abs(determinant)

    gram = gram_matrix(matrix)
    gram_determinant = bareiss_determinant(gram)
    if gram_determinant != determinant * determinant:
        raise VerificationError("det(A A^T) != det(A)^2")

    modular_checks: list[dict[str, int]] = []
    for prime in contract.modular_primes:
        observed = modular_determinant(matrix, prime)
        expected = determinant % prime
        if observed != expected:
            raise VerificationError(
                f"modular determinant disagrees for prime {prime}"
            )
        observed_gram = modular_determinant(gram, prime)
        expected_gram = determinant * determinant % prime
        if observed_gram != expected_gram:
            raise VerificationError(
                f"modular Gram determinant disagrees for prime {prime}"
            )
        modular_checks.append(
            {
                "prime": prime,
                "determinant_mod_prime": observed,
                "gram_determinant_mod_prime": observed_gram,
            }
        )

    bound_squared = contract.order ** contract.order
    determinant_squared = determinant * determinant
    if determinant_squared > bound_squared:
        raise VerificationError("Hadamard bound violated")

    ehlich_bound_squared = contract.order_specific_bound_squared
    if determinant_squared > ehlich_bound_squared:
        raise VerificationError("order-23 Ehlich bound violated")

    barba_bound_squared = (
        (2 * contract.order - 1)
        * (contract.order - 1) ** (contract.order - 1)
    )
    if determinant_squared > barba_bound_squared:
        raise VerificationError("Barba bound for odd order violated")

    required_factor = 1 << (contract.order - 1)
    if absolute % required_factor:
        raise VerificationError(
            f"determinant is not divisible by 2^{contract.order - 1}"
        )

    modulus_product = 1
    for prime in contract.modular_primes:
        modulus_product *= prime
    if modulus_product * modulus_product <= 4 * bound_squared:
        raise VerificationError(
            "modular witnesses do not uniquely determine a bounded determinant"
        )

    normalized = normalize_signs(matrix)
    normalized_bytes = matrix_text(normalized).encode("ascii")
    hadamard_ratio_squared_ppb = (
        determinant_squared * 1_000_000_000 // bound_squared
    )
    barba_ratio_squared_ppb = (
        determinant_squared * 1_000_000_000 // barba_bound_squared
    )
    ehlich_ratio_squared_ppb = (
        determinant_squared * 1_000_000_000 // ehlich_bound_squared
    )

    receipt: dict[str, Any] = {
        "receipt_schema_version": 2,
        "challenge_id": contract.challenge_id,
        "contract_sha256": contract.sha256,
        "verifier_version": VERIFIER_VERSION,
        "matrix": {
            "order": contract.order,
            "raw_sha256": hashlib.sha256(raw_bytes).hexdigest(),
            "sign_normalized_sha256": hashlib.sha256(normalized_bytes).hexdigest(),
            "normalization_scope": "row-and-column-signs-only",
        },
        "score": {
            "ranking_quantity": "absolute_determinant",
            "determinant": str(determinant),
            "absolute_determinant": str(absolute),
            "determinant_squared": str(determinant_squared),
            "hadamard_bound_squared": str(bound_squared),
            "hadamard_ratio_squared_parts_per_billion": (
                hadamard_ratio_squared_ppb
            ),
            "ehlich_bound_squared": str(ehlich_bound_squared),
            "ehlich_ratio_squared_parts_per_billion": ehlich_ratio_squared_ppb,
            "barba_bound_squared": str(barba_bound_squared),
            "barba_ratio_squared_parts_per_billion": barba_ratio_squared_ppb,
        },
        "checks": {
            "entry_domain": "passed",
            "exact_bareiss": "passed",
            "gram_identity": "passed",
            "gram_determinant": str(gram_determinant),
            "modular_determinants": modular_checks,
            "modular_modulus_product": str(modulus_product),
            "crt_unique_reconstruction": "passed",
            "order_23_ehlich_bound": "passed",
            "barba_bound": "passed",
            "hadamard_bound": "passed",
            "power_of_two_divisibility": "passed",
        },
        "claim_boundary": (
            "This receipt verifies the matrix and its exact score only. "
            "It does not establish global optimality or a world record."
        ),
    }

    receipt_bytes_without_hash = canonical_json_bytes(receipt)
    receipt["receipt_sha256"] = hashlib.sha256(
        receipt_bytes_without_hash
    ).hexdigest()
    return VerifiedMatrix(matrix=matrix, receipt=receipt)

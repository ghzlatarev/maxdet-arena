#!/usr/bin/env python3
"""Deterministic polar/isospectral Douglas--Rachford portal pilot.

The discrete set is C = {-1,+1}^{23x23}.  For a prescribed, descending
singular-value vector ``s``, the continuous set is

    D_s = {U diag(s) V^T : U,V orthogonal}.

Projection onto C is entrywise sign (zero maps to +1).  Projection onto D_s is
``U diag(s) V^T`` for an SVD ``X = U diag(_) V^T``.  The latter is a nearest
Frobenius-norm projection by the von Neumann trace inequality.

Two Douglas--Rachford orders are sampled:

    sign first:      X <- X + lambda(P_D(2 P_C(X) - X) - P_C(X))
    spectrum first:  X <- X + lambda(P_C(2 P_D(X) - X) - P_D(X))

Every distinct sign shadow that the pilot observes is scored with the trusted
Python Bareiss implementation.  Floating point only guides the search.  Ties
and wins are subsequently passed through ``./arena verify``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from maxdet.contract import load_contract, load_matrix
from maxdet.exact import bareiss_determinant, matrix_text
from maxdet.frontier import effective_frontier

ORDER = 23
ENTRIES = ORDER * ORDER
DEFAULT_AUDIT = (
    REPOSITORY_ROOT
    / "runs/direct-search/frontier-factor-class-expansion-20260728"
    / "final-h-equivalence-audit.json"
)
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT / "runs/direct-search/polar-dr-portals-20260729"
)
DEFAULT_CLASSIFIER = Path("/tmp/maxdet-h-audit-20260729/bin/python")
Matrix = tuple[tuple[int, ...], ...]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repository_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def immutable_matrix(matrix: Sequence[Sequence[int]]) -> Matrix:
    result = tuple(tuple(int(value) for value in row) for row in matrix)
    if len(result) != ORDER or any(len(row) != ORDER for row in result):
        raise ValueError("matrix is not 23 by 23")
    if any(value not in (-1, 1) for row in result for value in row):
        raise ValueError("matrix has an entry outside {-1,+1}")
    return result


def matrix_bits(matrix: Matrix) -> int:
    key = 0
    for row in matrix:
        for value in row:
            key = (key << 1) | (value > 0)
    return key


def matrix_hash(matrix: Matrix) -> str:
    return sha256_bytes(matrix_text(matrix).encode("ascii"))


def hamming_distance(first: Matrix, second: Matrix) -> int:
    return sum(
        left != right
        for left_row, right_row in zip(first, second)
        for left, right in zip(left_row, right_row)
    )


def project_sign(array: np.ndarray) -> Matrix:
    """Nearest sign matrix, with the deterministic convention sign(0)=+1."""
    if array.shape != (ORDER, ORDER):
        raise ValueError("projection input is not 23 by 23")
    signs = np.where(array < 0.0, -1, 1)
    return immutable_matrix(signs.tolist())


def project_spectrum(array: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Nearest matrix with the prescribed singular values."""
    if array.shape != (ORDER, ORDER) or target.shape != (ORDER,):
        raise ValueError("bad isospectral projection shape")
    if not np.all(np.isfinite(array)):
        raise FloatingPointError("non-finite Douglas--Rachford state")
    if np.any(target < 0.0) or np.any(target[:-1] < target[1:] - 1e-12):
        raise ValueError("target singular values must be nonnegative and descending")
    left, _, right_transpose = np.linalg.svd(array, full_matrices=False)
    return (left * target) @ right_transpose


def matrix_array(matrix: Matrix) -> np.ndarray:
    return np.asarray(matrix, dtype=np.float64)


def frontier_spectrum(matrix: Matrix) -> np.ndarray:
    return np.linalg.svd(matrix_array(matrix), compute_uv=False)


def flat_spectrum() -> np.ndarray:
    return np.full(ORDER, math.sqrt(ORDER), dtype=np.float64)


def ehlich_ideal_gram() -> np.ndarray:
    """Return G*=24I-J+4B for the disjoint K3 + 5 K4 partition."""
    gram = 24.0 * np.eye(ORDER) - np.ones((ORDER, ORDER))
    start = 0
    for size in (3, 4, 4, 4, 4, 4):
        stop = start + size
        gram[start:stop, start:stop] += 4.0 * (
            np.ones((size, size)) - np.eye(size)
        )
        start = stop
    if start != ORDER:
        raise AssertionError("Ehlich block partition does not have order 23")
    return gram


def ehlich_spectrum() -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh(ehlich_ideal_gram())
    if np.min(eigenvalues) <= 0.0:
        raise AssertionError("Ehlich ideal Gram is not positive definite")
    return np.sort(np.sqrt(eigenvalues))[::-1]


def normalized_homotopy(
    source: np.ndarray, destination: np.ndarray, fraction: float
) -> np.ndarray:
    """Interpolate spectra while preserving the sign-matrix Frobenius norm."""
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("homotopy fraction is outside [0,1]")
    target = (1.0 - fraction) * source + fraction * destination
    target *= ORDER / np.linalg.norm(target)
    return np.sort(target)[::-1]


def spectrum_for_iteration(
    mode: str,
    iteration: int,
    iterations: int,
    source: np.ndarray,
    flat: np.ndarray,
    ehlich: np.ndarray,
) -> np.ndarray:
    if mode == "frontier":
        return source
    if mode == "flat":
        return flat
    if mode == "ehlich":
        return ehlich
    denominator = max(1, int(math.ceil(0.8 * iterations)) - 1)
    fraction = min(1.0, iteration / denominator)
    if mode == "frontier-to-flat":
        return normalized_homotopy(source, flat, fraction)
    if mode == "frontier-to-ehlich":
        return normalized_homotopy(source, ehlich, fraction)
    raise ValueError(f"unknown spectrum mode: {mode}")


def derived_seed(base: int, *parts: object) -> int:
    material = "|".join([str(base), *(str(part) for part in parts)])
    return int.from_bytes(
        hashlib.sha256(material.encode("utf-8")).digest()[:8], "little"
    )


def canonical_qr(array: np.ndarray) -> np.ndarray:
    q_matrix, r_matrix = np.linalg.qr(array)
    diagonal = np.diag(r_matrix)
    signs = np.where(diagonal < 0.0, -1.0, 1.0)
    return q_matrix * signs


def perturb(
    source: Matrix, randomizer: np.random.Generator, restart: int
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return a deterministic continuous departure from one discrete portal."""
    base = matrix_array(source)
    family = restart % 3
    level = restart // 3
    if family == 0:
        sigma = 0.12 + 0.12 * min(level, 3)
        state = base + sigma * randomizer.standard_normal(base.shape)
        metadata = {"kind": "gaussian", "sigma": sigma}
    elif family == 1:
        strength = 0.08 + 0.08 * min(level, 3)
        left_noise = randomizer.standard_normal(base.shape)
        right_noise = randomizer.standard_normal(base.shape)
        left_skew = left_noise - left_noise.T
        right_skew = right_noise - right_noise.T
        scale = math.sqrt(ORDER)
        left = canonical_qr(np.eye(ORDER) + strength * left_skew / scale)
        right = canonical_qr(np.eye(ORDER) + strength * right_skew / scale)
        state = left @ base @ right.T
        state += 0.02 * randomizer.standard_normal(base.shape)
        metadata = {
            "kind": "isospectral_rotation",
            "strength": strength,
            "gaussian_sigma": 0.02,
        }
    else:
        kick_size = 4 + 4 * min(level, 5)
        kicked = base.copy()
        coordinates = randomizer.choice(ENTRIES, size=kick_size, replace=False)
        kicked.flat[coordinates] *= -1.0
        sigma = 0.08 + 0.04 * min(level, 3)
        state = kicked + sigma * randomizer.standard_normal(base.shape)
        metadata = {
            "kind": "entry_kick_gaussian",
            "entry_kick": kick_size,
            "sigma": sigma,
        }
    metadata["initial_shadow_hamming"] = hamming_distance(
        source, project_sign(state)
    )
    return state, metadata


@dataclass(frozen=True)
class Candidate:
    score: int
    key: int
    matrix: Matrix
    context: dict[str, Any]


class CandidateRegistry:
    """Deduplicate globally and exact-score every distinct observed matrix."""

    def __init__(self, frontier: int, per_class_limit: int) -> None:
        self.frontier = frontier
        self.per_class_limit = per_class_limit
        self.scores: dict[int, int] = {}
        self.observations = 0
        self.duplicates = 0
        self.distinct = 0
        self.exact_bareiss_scores = 0
        self.raw_observations = 0
        self.line_observations = 0
        self.raw_distinct = 0
        self.line_distinct = 0
        self.special: dict[int, Candidate] = {}
        self.class_top: dict[str, list[Candidate]] = {}
        self.best_subfrontier: Candidate | None = None
        self.best_overall: Candidate | None = None

    def observe(
        self, matrix: Matrix, context: dict[str, Any], category: str = "raw"
    ) -> tuple[int, bool]:
        self.observations += 1
        if category == "raw":
            self.raw_observations += 1
        elif category == "line":
            self.line_observations += 1
        else:
            raise ValueError(f"unknown candidate category: {category}")
        key = matrix_bits(matrix)
        existing = self.scores.get(key)
        if existing is not None:
            self.duplicates += 1
            return existing, False
        score = abs(bareiss_determinant(matrix))
        self.scores[key] = score
        self.distinct += 1
        self.exact_bareiss_scores += 1
        if category == "raw":
            self.raw_distinct += 1
        else:
            self.line_distinct += 1
        candidate = Candidate(score, key, matrix, dict(context))
        if (
            self.best_overall is None
            or (candidate.score, -candidate.key)
            > (self.best_overall.score, -self.best_overall.key)
        ):
            self.best_overall = candidate
        if score >= self.frontier:
            self.special[key] = candidate
        else:
            if (
                self.best_subfrontier is None
                or (candidate.score, -candidate.key)
                > (self.best_subfrontier.score, -self.best_subfrontier.key)
            ):
                self.best_subfrontier = candidate
            certificate = str(context["ht_certificate"])
            if category == "raw":
                retained = self.class_top.setdefault(certificate, [])
                retained.append(candidate)
                retained.sort(key=lambda item: (-item.score, item.key))
                del retained[self.per_class_limit :]
        return score, True


def exact_line_ascent(
    start: Candidate,
    registry: CandidateRegistry,
    max_sweeps: int,
    run_index: int,
) -> tuple[Candidate, dict[str, Any]]:
    """Exact best-row/column ascent using exact rational inverse signs."""
    import sympy

    current = start
    accepted = 0
    proposals = 0
    for sweep in range(max_sweeps):
        inverse = sympy.Matrix(current.matrix).inv(method="DM")
        best = current
        for row in range(ORDER):
            candidate_rows = [list(values) for values in current.matrix]
            for column in range(ORDER):
                value = inverse[column, row]
                if value > 0:
                    candidate_rows[row][column] = 1
                elif value < 0:
                    candidate_rows[row][column] = -1
            matrix = immutable_matrix(candidate_rows)
            context = {
                **start.context,
                "stage": "exact_line_ascent",
                "line_kind": "row",
                "line_index": row,
                "line_sweep": sweep,
                "line_run_index": run_index,
            }
            score, _ = registry.observe(matrix, context, "line")
            proposals += 1
            proposed = Candidate(score, matrix_bits(matrix), matrix, context)
            if (score, -proposed.key) > (best.score, -best.key):
                best = proposed
        for column in range(ORDER):
            candidate_rows = [list(values) for values in current.matrix]
            for row in range(ORDER):
                value = inverse[column, row]
                if value > 0:
                    candidate_rows[row][column] = 1
                elif value < 0:
                    candidate_rows[row][column] = -1
            matrix = immutable_matrix(candidate_rows)
            context = {
                **start.context,
                "stage": "exact_line_ascent",
                "line_kind": "column",
                "line_index": column,
                "line_sweep": sweep,
                "line_run_index": run_index,
            }
            score, _ = registry.observe(matrix, context, "line")
            proposals += 1
            proposed = Candidate(score, matrix_bits(matrix), matrix, context)
            if (score, -proposed.key) > (best.score, -best.key):
                best = proposed
        if best.score <= current.score:
            break
        current = best
        accepted += 1
    return current, {
        "start_score": str(start.score),
        "final_score": str(current.score),
        "strict_line_moves": accepted,
        "line_proposals": proposals,
        "start_matrix_sha256": matrix_hash(start.matrix),
        "final_matrix_sha256": matrix_hash(current.matrix),
    }


def load_portals(audit_path: Path, frontier: int) -> list[dict[str, Any]]:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("expected_absolute_determinant") != frontier:
        raise RuntimeError("frontier audit has the wrong expected score")
    classes = audit.get("ht_classes")
    if not isinstance(classes, list) or len(classes) != 8:
        raise RuntimeError("frontier audit does not contain exactly 8 H/HT classes")
    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    portals: list[dict[str, Any]] = []
    for class_index, item in enumerate(classes):
        members = item.get("members")
        if not isinstance(members, list) or not members:
            raise RuntimeError("H/HT class has no member")
        path = REPOSITORY_ROOT / members[0]["path"]
        raw, loaded = load_matrix(path, contract)
        matrix = immutable_matrix(loaded)
        score = abs(bareiss_determinant(matrix))
        if score != frontier:
            raise RuntimeError(f"{path}: portal score {score} != frontier {frontier}")
        portals.append(
            {
                "class_index": class_index,
                "ht_certificate": item["certificate_sha256"],
                "path": path,
                "path_label": repository_path(path),
                "raw_sha256": sha256_bytes(raw),
                "matrix_sha256": matrix_hash(matrix),
                "score": score,
                "matrix": matrix,
            }
        )
    return portals


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def write_matrix(path: Path, matrix: Matrix) -> None:
    atomic_write_text(path, matrix_text(matrix))


def append_event(path: Path, event: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as output:
        output.write(json.dumps(event, sort_keys=True) + "\n")
        output.flush()


def run_trajectory(
    portal: dict[str, Any],
    mode: str,
    order: str,
    restart: int,
    arguments: argparse.Namespace,
    registry: CandidateRegistry,
    event_path: Path,
    campaign_started: float,
    next_heartbeat: list[float],
) -> dict[str, Any]:
    seed = derived_seed(
        arguments.seed, portal["ht_certificate"], mode, order, restart
    )
    randomizer = np.random.default_rng(seed)
    state, perturbation = perturb(portal["matrix"], randomizer, restart)
    source_spectrum = frontier_spectrum(portal["matrix"])
    flat = flat_spectrum()
    ehlich = ehlich_spectrum()
    relaxation = 0.8 if restart % 2 == 0 else 0.9
    local_seen: set[int] = set()
    local_distinct = 0
    local_duplicates = 0
    best_score = -1
    best_hash = ""
    last_local_new = 0
    escape_kicks = 0
    numerical_resets = 0

    def observe(matrix: Matrix, iteration: int, shadow: str) -> None:
        nonlocal local_distinct, local_duplicates, best_score, best_hash
        nonlocal last_local_new
        key = matrix_bits(matrix)
        if key in local_seen:
            local_duplicates += 1
        else:
            local_seen.add(key)
            local_distinct += 1
            last_local_new = iteration
        context = {
            "stage": "douglas_rachford",
            "class_index": portal["class_index"],
            "ht_certificate": portal["ht_certificate"],
            "source_path": portal["path_label"],
            "spectrum_mode": mode,
            "dr_order": order,
            "restart": restart,
            "trajectory_seed": seed,
            "iteration": iteration,
            "shadow": shadow,
        }
        score, _ = registry.observe(matrix, context, "raw")
        candidate_hash = matrix_hash(matrix)
        if (score, candidate_hash) > (best_score, best_hash):
            best_score = score
            best_hash = candidate_hash

    for iteration in range(arguments.iterations):
        target = spectrum_for_iteration(
            mode,
            iteration,
            arguments.iterations,
            source_spectrum,
            flat,
            ehlich,
        )
        try:
            if order == "sign-first":
                sign = project_sign(state)
                observe(sign, iteration, "P_C(X)")
                spectral = project_spectrum(
                    2.0 * matrix_array(sign) - state, target
                )
                observe(project_sign(spectral), iteration, "P_C(P_D(R_C(X)))")
                state = state + relaxation * (spectral - matrix_array(sign))
                observe(project_sign(state), iteration, "P_C(X_next)")
            elif order == "spectrum-first":
                spectral = project_spectrum(state, target)
                observe(project_sign(spectral), iteration, "P_C(P_D(X))")
                sign = project_sign(2.0 * spectral - state)
                observe(sign, iteration, "P_C(R_D(X))")
                state = state + relaxation * (matrix_array(sign) - spectral)
                observe(project_sign(state), iteration, "P_C(X_next)")
            else:
                raise ValueError(f"unknown Douglas--Rachford order: {order}")
        except (FloatingPointError, np.linalg.LinAlgError):
            numerical_resets += 1
            state, _ = perturb(portal["matrix"], randomizer, restart + 97)
            continue

        if iteration - last_local_new >= arguments.stagnation_iterations:
            kick_scale = 0.04 + 0.02 * (escape_kicks % 3)
            state += kick_scale * randomizer.standard_normal(state.shape)
            escape_kicks += 1
            last_local_new = iteration

        now = time.monotonic()
        if arguments.heartbeat_seconds and now >= next_heartbeat[0]:
            append_event(
                event_path,
                {
                    "event": "heartbeat",
                    "elapsed_seconds": round(now - campaign_started, 3),
                    "distinct_candidates": registry.distinct,
                    "duplicate_observations": registry.duplicates,
                    "exact_bareiss_scores": registry.exact_bareiss_scores,
                    "trajectory": {
                        "class_index": portal["class_index"],
                        "spectrum_mode": mode,
                        "dr_order": order,
                        "restart": restart,
                        "iteration": iteration,
                    },
                },
            )
            next_heartbeat[0] = now + arguments.heartbeat_seconds

    final_target = spectrum_for_iteration(
        mode,
        arguments.iterations - 1,
        arguments.iterations,
        source_spectrum,
        flat,
        ehlich,
    )
    final_spectral = project_spectrum(state, final_target)
    return {
        "class_index": portal["class_index"],
        "ht_certificate": portal["ht_certificate"],
        "source_path": portal["path_label"],
        "spectrum_mode": mode,
        "dr_order": order,
        "restart": restart,
        "trajectory_seed": seed,
        "relaxation": relaxation,
        "iterations": arguments.iterations,
        "shadow_observations": arguments.iterations * 3,
        "local_distinct_shadows": local_distinct,
        "local_duplicate_shadows": local_duplicates,
        "best_observed_score": str(best_score),
        "best_observed_matrix_sha256": best_hash,
        "escape_kicks": escape_kicks,
        "numerical_resets": numerical_resets,
        "perturbation": perturbation,
        "final_spectral_projection_residual_frobenius": float(
            np.linalg.norm(
                np.linalg.svd(final_spectral, compute_uv=False) - final_target
            )
        ),
    }


def dependency_versions() -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for package in ("numpy", "scipy", "sympy", "pynauty"):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = None
    return result


def projection_self_tests() -> dict[str, Any]:
    sign_input = np.ones((ORDER, ORDER))
    sign_input[0, 0] = -0.1
    sign_input[0, 1] = 0.0
    sign = project_sign(sign_input)
    if sign[0][0] != -1 or sign[0][1] != 1:
        raise AssertionError("deterministic sign projection failed")

    randomizer = np.random.default_rng(230729)
    array = randomizer.standard_normal((ORDER, ORDER))
    target = ehlich_spectrum()
    projected = project_spectrum(array, target)
    observed = np.linalg.svd(projected, compute_uv=False)
    spectral_error = float(np.max(np.abs(observed - target)))
    if spectral_error > 2e-11:
        raise AssertionError("isospectral projection self-test failed")

    source_path = REPOSITORY_ROOT / "references/orrick-et-al-2003/matrix.txt"
    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    _, source_rows = load_matrix(source_path, contract)
    source = immutable_matrix(source_rows)
    expected = 2_779_447_296_000_000
    if abs(bareiss_determinant(source)) != expected:
        raise AssertionError("Bareiss reference-score self-test failed")

    def deterministic_trace(seed: int) -> list[int]:
        rng = np.random.default_rng(seed)
        state, _ = perturb(source, rng, 1)
        result: list[int] = []
        target_values = flat_spectrum()
        for _ in range(6):
            discrete = project_sign(state)
            result.append(matrix_bits(discrete))
            spectral = project_spectrum(
                2.0 * matrix_array(discrete) - state, target_values
            )
            state += 0.8 * (spectral - matrix_array(discrete))
        return result

    if deterministic_trace(91) != deterministic_trace(91):
        raise AssertionError("deterministic trajectory self-test failed")
    frontier_values = frontier_spectrum(source)
    if not np.allclose(
        normalized_homotopy(frontier_values, flat_spectrum(), 0.0),
        frontier_values,
        rtol=0.0,
        atol=1e-12,
    ):
        raise AssertionError("homotopy start self-test failed")
    if not np.allclose(
        normalized_homotopy(frontier_values, flat_spectrum(), 1.0),
        flat_spectrum(),
        rtol=0.0,
        atol=1e-12,
    ):
        raise AssertionError("homotopy end self-test failed")
    ideal = ehlich_ideal_gram()
    if not np.all(np.diag(ideal) == 23.0):
        raise AssertionError("Ehlich target diagonal self-test failed")
    ideal_off_diagonal = ideal[~np.eye(ORDER, dtype=bool)]
    if set(np.unique(ideal_off_diagonal)) != {-1.0, 3.0}:
        raise AssertionError("Ehlich target off-diagonal self-test failed")
    return {
        "passed": True,
        "sign_zero_convention": "+1",
        "isospectral_max_singular_value_error": spectral_error,
        "reference_absolute_determinant": str(expected),
        "deterministic_trace_length": 6,
        "ehlich_ideal_gram_determinant": str(
            bareiss_determinant(
                [[int(value) for value in row] for row in ideal.tolist()]
            )
        ),
    }


def verify_special_candidates(
    registry: CandidateRegistry, output_dir: Path
) -> list[dict[str, Any]]:
    verified: list[dict[str, Any]] = []
    for candidate in sorted(
        registry.special.values(), key=lambda item: (-item.score, item.key)
    ):
        category = "wins" if candidate.score > registry.frontier else "frontier-ties"
        digest = matrix_hash(candidate.matrix)
        matrix_path = output_dir / category / f"{digest}.matrix.txt"
        receipt_path = output_dir / category / f"{digest}.receipt.json"
        write_matrix(matrix_path, candidate.matrix)
        completed = subprocess.run(
            [
                str(REPOSITORY_ROOT / "arena"),
                "verify",
                str(matrix_path),
                "--json",
                str(receipt_path),
            ],
            cwd=REPOSITORY_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"arena verification failed for {matrix_path}: "
                f"{completed.stderr.strip()}"
            )
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        observed = int(receipt["score"]["absolute_determinant"])
        if observed != candidate.score:
            raise RuntimeError("arena receipt disagrees with Bareiss candidate score")
        verified.append(
            {
                "category": category,
                "absolute_determinant": str(candidate.score),
                "matrix_path": repository_path(matrix_path),
                "matrix_sha256": digest,
                "receipt_path": repository_path(receipt_path),
                "receipt_sha256": receipt["receipt_sha256"],
                "receipt_file_sha256": sha256_file(receipt_path),
                "normalized_sha256": receipt["matrix"]["sign_normalized_sha256"],
                "first_context": candidate.context,
            }
        )
    return verified


def classify_frontier_ties(
    verified: Sequence[dict[str, Any]],
    classifier_python: Path | None,
    frontier: int,
    output_dir: Path,
) -> dict[str, Any]:
    ties = [
        REPOSITORY_ROOT / item["matrix_path"]
        for item in verified
        if item["category"] == "frontier-ties"
    ]
    if not ties:
        return {"status": "not-needed", "tie_count": 0}
    if classifier_python is None or not classifier_python.exists():
        return {
            "status": "unavailable",
            "tie_count": len(ties),
            "classifier_python": (
                None if classifier_python is None else str(classifier_python)
            ),
        }
    audit_path = output_dir / "frontier-tie-h-equivalence.json"
    completed = subprocess.run(
        [
            str(classifier_python),
            str(REPOSITORY_ROOT / "research/h_equivalence_audit.py"),
            "--expected-absolute-determinant",
            str(frontier),
            "--output",
            str(audit_path),
            *(str(path) for path in ties),
        ],
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"H/HT classifier failed: {completed.stderr.strip()}"
        )
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    return {
        "status": "complete",
        "tie_count": len(ties),
        "classifier_python": str(classifier_python),
        "classifier_python_sha256": sha256_file(classifier_python),
        "pynauty_version": audit.get("pynauty_version"),
        "h_class_count": audit["h_class_count"],
        "ht_class_count": audit["ht_class_count"],
        "gram_class_count": audit["gram_class_count"],
        "audit_path": repository_path(audit_path),
        "audit_sha256": sha256_file(audit_path),
    }


def candidate_record(candidate: Candidate | None, path: Path | None) -> Any:
    if candidate is None:
        return None
    result = {
        "absolute_determinant": str(candidate.score),
        "matrix_sha256": matrix_hash(candidate.matrix),
        "first_context": candidate.context,
    }
    if path is not None:
        result["path"] = repository_path(path)
        result["raw_sha256"] = sha256_file(path)
    return result


def git_metadata() -> dict[str, Any]:
    def command(*arguments: str) -> str:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=REPOSITORY_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return completed.stdout.strip()

    return {
        "commit": command("rev-parse", "HEAD"),
        "branch": command("branch", "--show-current"),
        "status_short": command("status", "--short"),
    }


def run_campaign(arguments: argparse.Namespace) -> int:
    output_dir = arguments.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"refusing nonempty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    event_path = output_dir / "events.jsonl"
    started_wall = utc_now()
    started = time.monotonic()
    self_tests = projection_self_tests()
    contract = load_contract(REPOSITORY_ROOT / "challenge.json")
    frontier_record = effective_frontier(REPOSITORY_ROOT, contract)
    frontier = frontier_record.absolute_determinant
    portals = load_portals(arguments.audit.resolve(), frontier)
    registry = CandidateRegistry(
        frontier,
        per_class_limit=max(arguments.postprocess_per_class * 3, 6),
    )
    append_event(
        event_path,
        {
            "event": "start",
            "started_at": started_wall,
            "seed": arguments.seed,
            "portal_count": len(portals),
            "frontier": str(frontier),
            "self_tests": self_tests,
        },
    )

    # Register known portals for accounting and verification, but report them as
    # inputs rather than as search discoveries.
    for portal in portals:
        registry.observe(
            portal["matrix"],
            {
                "stage": "input_portal",
                "class_index": portal["class_index"],
                "ht_certificate": portal["ht_certificate"],
                "source_path": portal["path_label"],
            },
            "raw",
        )

    trajectories: list[dict[str, Any]] = []
    next_heartbeat = [started + arguments.heartbeat_seconds]
    for portal in portals:
        for mode in arguments.spectrum_modes:
            for order in arguments.dr_orders:
                for restart in range(arguments.restarts):
                    trajectory = run_trajectory(
                        portal,
                        mode,
                        order,
                        restart,
                        arguments,
                        registry,
                        event_path,
                        started,
                        next_heartbeat,
                    )
                    trajectories.append(trajectory)
                    append_event(
                        event_path,
                        {
                            "event": "trajectory_finished",
                            "elapsed_seconds": round(
                                time.monotonic() - started, 3
                            ),
                            **trajectory,
                        },
                    )

    postprocess_inputs: list[Candidate] = []
    for portal in portals:
        certificate = portal["ht_certificate"]
        postprocess_inputs.extend(
            registry.class_top.get(certificate, [])[
                : arguments.postprocess_per_class
            ]
        )
    unique_inputs: dict[int, Candidate] = {
        candidate.key: candidate for candidate in postprocess_inputs
    }
    line_runs: list[dict[str, Any]] = []
    line_outputs: dict[int, Candidate] = {}
    for run_index, candidate in enumerate(
        sorted(
            unique_inputs.values(), key=lambda item: (-item.score, item.key)
        )
    ):
        final, metadata = exact_line_ascent(
            candidate,
            registry,
            arguments.line_max_sweeps,
            run_index,
        )
        line_outputs[final.key] = final
        line_runs.append(metadata)

    best_path: Path | None = None
    if registry.best_subfrontier is not None:
        best_path = output_dir / "best-subfrontier.matrix.txt"
        write_matrix(best_path, registry.best_subfrontier.matrix)
    raw_top_paths: list[dict[str, Any]] = []
    raw_candidates = {
        item.key: item
        for retained in registry.class_top.values()
        for item in retained
    }
    for index, candidate in enumerate(
        sorted(raw_candidates.values(), key=lambda item: (-item.score, item.key))[
            : arguments.retain
        ]
    ):
        path = (
            output_dir
            / "top-subfrontier"
            / f"{index:03d}-{candidate.score}-{matrix_hash(candidate.matrix)}.matrix.txt"
        )
        write_matrix(path, candidate.matrix)
        raw_top_paths.append(candidate_record(candidate, path))
    line_output_paths: list[dict[str, Any]] = []
    for index, candidate in enumerate(
        sorted(line_outputs.values(), key=lambda item: (-item.score, item.key))[
            : arguments.retain
        ]
    ):
        path = (
            output_dir
            / "line-ascended"
            / f"{index:03d}-{candidate.score}-{matrix_hash(candidate.matrix)}.matrix.txt"
        )
        write_matrix(path, candidate.matrix)
        line_output_paths.append(candidate_record(candidate, path))

    verified = verify_special_candidates(registry, output_dir)
    classifier = classify_frontier_ties(
        verified,
        arguments.classifier_python,
        frontier,
        output_dir,
    )
    finished = time.monotonic()
    special_search_discoveries = [
        item
        for item in verified
        if item["first_context"].get("stage") != "input_portal"
    ]
    report = {
        "schema_version": 1,
        "engine": "polar-isospectral-douglas-rachford-portals-v1",
        "claim_boundary": (
            "Continuous projections guide search only. Every distinct sign "
            "shadow is Bareiss-scored; arena receipts verify all ties/wins. "
            "The Ehlich ideal spectrum is guidance, not a factorization claim."
        ),
        "complete": True,
        "started_at": started_wall,
        "finished_at": utc_now(),
        "elapsed_seconds": round(finished - started, 6),
        "frontier": {
            "absolute_determinant": str(frontier),
            "source": frontier_record.source,
        },
        "configuration": {
            "seed": arguments.seed,
            "iterations": arguments.iterations,
            "restarts": arguments.restarts,
            "spectrum_modes": arguments.spectrum_modes,
            "dr_orders": arguments.dr_orders,
            "stagnation_iterations": arguments.stagnation_iterations,
            "postprocess_per_class": arguments.postprocess_per_class,
            "line_max_sweeps": arguments.line_max_sweeps,
            "retain": arguments.retain,
        },
        "projection_definition": {
            "sign": "entrywise sign with sign(0)=+1",
            "isospectral": "if X=U diag(q) V^T, P_D(X)=U diag(target) V^T",
            "sign_first_update": (
                "X += lambda * (P_D(2*P_C(X)-X) - P_C(X))"
            ),
            "spectrum_first_update": (
                "X += lambda * (P_C(2*P_D(X)-X) - P_D(X))"
            ),
            "ehlich_ideal_gram": "24I-J+4B; B adjacency of disjoint K3 + 5 K4",
        },
        "self_tests": self_tests,
        "dependency_audit": dependency_versions(),
        "portals": [
            {
                key: value
                for key, value in portal.items()
                if key not in {"matrix", "path"}
            }
            for portal in portals
        ],
        "counts": {
            "portal_classes": len(portals),
            "trajectories": len(trajectories),
            "candidate_observations": registry.observations,
            "duplicate_observations": registry.duplicates,
            "distinct_candidates": registry.distinct,
            "exact_bareiss_scores": registry.exact_bareiss_scores,
            "raw_shadow_observations_including_inputs": registry.raw_observations,
            "raw_distinct_candidates_including_inputs": registry.raw_distinct,
            "line_candidate_observations": registry.line_observations,
            "line_distinct_candidates": registry.line_distinct,
            "input_portals": len(portals),
            "arena_verified_ties_and_wins": len(verified),
            "search_discovered_ties_and_wins": len(special_search_discoveries),
            "strict_wins": sum(
                item["category"] == "wins" for item in verified
            ),
            "postprocess_unique_inputs": len(unique_inputs),
            "postprocess_unique_outputs": len(line_outputs),
        },
        "best_overall_including_input_portals": candidate_record(
            registry.best_overall, None
        ),
        "best_subfrontier": candidate_record(
            registry.best_subfrontier, best_path
        ),
        "search_discovered_ties_and_wins": special_search_discoveries,
        "arena_verified_ties_and_wins": verified,
        "frontier_tie_classification": classifier,
        "top_subfrontier": raw_top_paths,
        "line_ascended": line_output_paths,
        "line_ascent_runs": line_runs,
        "trajectories": trajectories,
    }
    report_path = output_dir / "report.json"
    atomic_write_text(report_path, json.dumps(report, indent=2, sort_keys=True) + "\n")
    append_event(
        event_path,
        {
            "event": "finished",
            "elapsed_seconds": round(finished - started, 3),
            "report_sha256": sha256_file(report_path),
            "best_subfrontier": (
                None
                if registry.best_subfrontier is None
                else str(registry.best_subfrontier.score)
            ),
            "strict_wins": report["counts"]["strict_wins"],
        },
    )

    provenance = {
        "schema_version": 1,
        "report_path": repository_path(report_path),
        "report_sha256": sha256_file(report_path),
        "tool_path": repository_path(Path(__file__)),
        "tool_sha256": sha256_file(Path(__file__)),
        "challenge_path": "challenge.json",
        "challenge_sha256": sha256_file(REPOSITORY_ROOT / "challenge.json"),
        "input_audit_path": repository_path(arguments.audit),
        "input_audit_sha256": sha256_file(arguments.audit),
        "event_log_path": repository_path(event_path),
        "event_log_sha256": sha256_file(event_path),
        "command": [sys.executable, *sys.argv],
        "python": {
            "executable": sys.executable,
            "version": sys.version,
            "platform": platform.platform(),
        },
        "dependencies": dependency_versions(),
        "git": git_metadata(),
    }
    provenance_path = output_dir / "provenance.json"
    atomic_write_text(
        provenance_path, json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    for path in (report_path, provenance_path):
        atomic_write_text(
            path.with_suffix(path.suffix + ".sha256"),
            f"{sha256_file(path)}  {path.name}\n",
        )
    print(json.dumps(
        {
            "report": repository_path(report_path),
            "distinct_candidates": registry.distinct,
            "duplicate_observations": registry.duplicates,
            "best_subfrontier": (
                None
                if registry.best_subfrontier is None
                else str(registry.best_subfrontier.score)
            ),
            "search_discovered_ties_and_wins": len(special_search_discoveries),
            "strict_wins": report["counts"]["strict_wins"],
        },
        sort_keys=True,
    ))
    return 0


def positive_integer(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def nonnegative_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed < 0.0:
        raise argparse.ArgumentTypeError("value must be finite and nonnegative")
    return parsed


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=23_072_923)
    parser.add_argument("--iterations", type=positive_integer, default=320)
    parser.add_argument("--restarts", type=positive_integer, default=6)
    parser.add_argument(
        "--spectrum-modes",
        nargs="+",
        choices=(
            "frontier",
            "flat",
            "ehlich",
            "frontier-to-flat",
            "frontier-to-ehlich",
        ),
        default=(
            "frontier",
            "flat",
            "ehlich",
            "frontier-to-flat",
            "frontier-to-ehlich",
        ),
    )
    parser.add_argument(
        "--dr-orders",
        nargs="+",
        choices=("sign-first", "spectrum-first"),
        default=("sign-first", "spectrum-first"),
    )
    parser.add_argument("--stagnation-iterations", type=positive_integer, default=48)
    parser.add_argument("--postprocess-per-class", type=positive_integer, default=3)
    parser.add_argument("--line-max-sweeps", type=positive_integer, default=12)
    parser.add_argument("--retain", type=positive_integer, default=48)
    parser.add_argument("--heartbeat-seconds", type=nonnegative_float, default=30.0)
    parser.add_argument(
        "--classifier-python",
        type=Path,
        default=DEFAULT_CLASSIFIER if DEFAULT_CLASSIFIER.exists() else None,
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run projection/determinism tests and exit without creating output",
    )
    arguments = parser.parse_args()
    if arguments.self_test:
        print(json.dumps(projection_self_tests(), indent=2, sort_keys=True))
        raise SystemExit(0)
    return arguments


if __name__ == "__main__":
    raise SystemExit(run_campaign(parse_arguments()))

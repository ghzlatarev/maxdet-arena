"""Small strict-JSON helpers for trusted arena inputs."""

from __future__ import annotations

import json
from typing import Any


class StrictJsonError(ValueError):
    """JSON is malformed, ambiguous, or not encoded as UTF-8."""


def loads_strict_json(raw: bytes) -> Any:
    """Decode UTF-8 JSON while rejecting duplicates and non-standard numbers."""

    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise StrictJsonError(f"not valid UTF-8: {exc}") from exc

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise StrictJsonError(f"duplicate object key: {key}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise StrictJsonError(f"non-standard numeric constant: {value}")

    try:
        return json.loads(
            text,
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except StrictJsonError:
        raise
    except (json.JSONDecodeError, RecursionError, ValueError) as exc:
        raise StrictJsonError(str(exc)) from exc

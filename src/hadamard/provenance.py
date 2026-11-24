"""Append-only provenance, branch artifacts, and deterministic JSON output."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


@dataclass(slots=True)
class Event:
    event: str
    branch_id: str | None = None
    payload: dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=utc_now)


class JsonlRecorder:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _default(value: Any) -> Any:
        if hasattr(value, "to_json"):
            return value.to_json()
        return repr(value)

    def write(self, event: Event | str, *, branch_id: str | None = None, **payload: Any) -> None:
        record = event if isinstance(event, Event) else Event(event, branch_id, dict(payload))
        with self.path.open("a", encoding="utf-8") as handle:
            json.dump(asdict(record), handle, sort_keys=True, default=self._default)
            handle.write("\n")


def write_json(path: str | Path, value: Mapping[str, Any] | list[Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, default=JsonlRecorder._default)
        handle.write("\n")

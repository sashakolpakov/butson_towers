#!/usr/bin/env python3
"""Run an isolated research command under a process-tree RSS/time watchdog."""

from __future__ import annotations

import argparse
import json
import math
import os
import signal
import subprocess
import time
from contextlib import suppress
from pathlib import Path

from hadamard.distributed import write_checkpoint
from hadamard.torsion.frontier import _process_tree_rss_records


def _process_group_exists(group_id: int) -> bool:
    try:
        os.killpg(group_id, 0)
    except ProcessLookupError:
        return False
    return True


def _cleanup_process_group(process: subprocess.Popen, grace_seconds: float = 3) -> None:
    """Stop the entire isolated group, including children of an exited root.

    Waiting only for the root process misses descendants that ignore TERM
    or outlive a successful parent.  Group membership remains available
    after the parent exits, so the grace period and final KILL apply to
    that group independently of the parent's return code.
    """

    group_id = process.pid
    with suppress(ProcessLookupError):
        os.killpg(group_id, signal.SIGTERM)
    deadline = time.monotonic() + grace_seconds
    while _process_group_exists(group_id) and time.monotonic() < deadline:
        process.poll()
        time.sleep(0.05)
    if _process_group_exists(group_id):
        with suppress(ProcessLookupError):
            os.killpg(group_id, signal.SIGKILL)
    process.wait()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-mb", type=int, default=512)
    parser.add_argument("--seconds", type=float, default=600)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.memory_mb <= 0 or args.seconds <= 0 or not math.isfinite(args.seconds):
        parser.error("a command and positive finite budgets are required")
    started = time.monotonic()
    peak = 0
    reason = None
    process = subprocess.Popen(command, start_new_session=True)
    try:
        while process.poll() is None:
            records = _process_tree_rss_records(process.pid)
            if not records and process.poll() is None:
                reason = "process-tree memory accounting unavailable"
                break
            peak = max(peak, sum(records.values()))
            if peak > args.memory_mb * 1024**2:
                reason = "process-tree memory limit"
                break
            if time.monotonic() - started > args.seconds:
                reason = "wall-clock limit"
                break
            time.sleep(0.1)
    finally:
        _cleanup_process_group(process)
    report = {
        "schema": "hadamard.bounded-command-report.v1",
        "command": command,
        "exit_code": process.returncode,
        "stop_reason": reason,
        "peak_process_tree_rss_bytes": peak,
        "elapsed_seconds": time.monotonic() - started,
        "memory_limit_mb": args.memory_mb,
        "time_limit_seconds": args.seconds,
        "classification_complete": False,
    }
    write_checkpoint(args.report, report)
    print(json.dumps(report, sort_keys=True), flush=True)
    return 2 if reason else (process.returncode or 0)


if __name__ == "__main__":
    raise SystemExit(main())

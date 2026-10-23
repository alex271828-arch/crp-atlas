"""Durable journals and atomic snapshots. Writers hold a per-pass file lock."""

from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path


@contextmanager
def lock(path, blocking=False):
    with Path(path).open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        except BlockingIOError:
            raise RuntimeError(f"Another process holds {path}") from None
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def atomic_json(path, record):
    path = Path(path)
    temporary = path.with_suffix(f".{os.getpid()}.tmp")
    with temporary.open("w") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def rows(path, repair_tail=False):
    """Read a journal, optionally removing a crash-torn final record under lock.

    A checkpoint is written only after a complete newline-terminated record is
    fsynced. A trailing partial write can therefore be replayed from the saved
    agenda. Malformed complete records are errors and are never silently lost.
    """
    path = Path(path)
    if not path.exists():
        return []
    result = []
    with path.open("r+b" if repair_tail else "rb") as handle:
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            if not line.endswith(b"\n"):
                if not repair_tail:
                    raise ValueError(
                        f"Incomplete journal tail: {path}; resume this pass to repair it"
                    )
                handle.seek(offset)
                handle.truncate()
                handle.flush()
                os.fsync(handle.fileno())
                break
            if line.strip():
                result.append(json.loads(line))
    return result


def append(path, record):
    with Path(path).open("a") as handle:
        handle.write(json.dumps(record) + "\n")
        handle.flush()
        os.fsync(handle.fileno())

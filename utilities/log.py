
from __future__ import annotations

import atexit
import os
import shlex
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path


DETAIL  = 5
INFO    = 20
STEP    = 25
STAGE   = 27
WARNING = 30
ERROR   = 40

_LEVEL_NAMES = {
    DETAIL: 'DETAIL', INFO: 'INFO', STEP: 'STEP',
    STAGE: 'STAGE', WARNING: 'WARNING', ERROR: 'ERROR',
}


_terminal_level = STEP
_file_level     = DETAIL


def terminal_level() -> int:
    return _terminal_level


def enabled(level: int) -> bool:
    """True if `level` reaches either sink — for guarding expensive formatting."""
    return level >= _terminal_level or level >= _file_level


_term_fh = None      
_log_fh  = None
_log_path: Path | None = None
_logs_dir: Path | None = None
_pending: list = []
_progress_open = False


def _terminal():
    return _term_fh if _term_fh is not None else sys.__stderr__


def _isatty(fh) -> bool:
    try:
        return bool(fh.isatty())
    except (AttributeError, OSError, ValueError):
        return False


def _emit(level: int, msg) -> None:
    global _progress_open
    text = str(msg)

    if level >= _terminal_level:
        fh = _terminal()
        try:
            if _progress_open:
                fh.write('\n')
                _progress_open = False
            fh.write(text + '\n')
            fh.flush()
        except (OSError, ValueError):
            pass

    if level >= _file_level:
        if _log_fh is None:
            _pending.append((level, text))
        else:
            try:
                _log_fh.write(text + '\n')
            except (OSError, ValueError):
                pass


def stage(msg='') -> None:
    _emit(STAGE, msg)


def step(msg='') -> None:
    _emit(STEP, msg)


def info(msg='') -> None:
    _emit(INFO, msg)


def detail(msg='') -> None:
    _emit(DETAIL, msg)


def warn(msg='') -> None:
    _emit(WARNING, f"Warning: {msg}")


def error(msg='') -> None:
    _emit(ERROR, msg)


def block(level, text) -> None:
    for line in str(text).splitlines() or ['']:
        _emit(level, line)


# ---------------------------------------------------------------------------
# Transient progress
# ---------------------------------------------------------------------------

def progress(msg) -> None:
    global _progress_open
    if _terminal_level > STEP:
        return
    fh = _terminal()
    if not _isatty(fh):
        return
    try:
        fh.write('\r' + str(msg))
        fh.flush()
        _progress_open = True
    except (OSError, ValueError):
        pass


def progress_end() -> None:
    global _progress_open
    if not _progress_open:
        return
    try:
        fh = _terminal()
        fh.write('\n')
        fh.flush()
    except (OSError, ValueError):
        pass
    _progress_open = False


# ---------------------------------------------------------------------------
# Tee: write to two streams simultaneously
# ---------------------------------------------------------------------------

class _Tee:
    def __init__(self, primary, secondary):
        self._primary   = primary
        self._secondary = secondary

    def write(self, data):
        self._primary.write(data)
        self._secondary.write(data)

    def flush(self):
        self._primary.flush()
        self._secondary.flush()

    # Make subprocess and other code that calls fileno() happy
    def fileno(self):
        return self._primary.fileno()

    def isatty(self):
        return self._primary.isatty()


_tee_thread = None
_saved_stderr_fd = None


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

def setup_logging(log_path) -> Path:
    global _log_fh, _log_path, _logs_dir, _term_fh, _tee_thread, _saved_stderr_fd

    log_path = Path(log_path)
    _log_path = log_path
    _logs_dir = log_path.parent / 'logs'
    _logs_dir.mkdir(parents=True, exist_ok=True)

    _log_fh = open(log_path, 'a', buffering=1)  
    _log_fh.write(
        f"\n# {'═' * 70}\n"
        f"# HiRepEat run log\n"
        f"# started : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"# command : {shlex.join(sys.argv)}\n"
        f"# log     : {log_path}\n"
        f"# tools   : {_logs_dir / 'tools.log'}\n"
        f"# {'═' * 70}\n"
    )

    _saved_stderr_fd = os.dup(2)
    r_fd, w_fd = os.pipe()
    os.dup2(w_fd, 2)
    os.close(w_fd)

    original_stderr = os.fdopen(_saved_stderr_fd, 'w', buffering=1, closefd=False)

    def _tee_fd2():
        while True:
            try:
                data = os.read(r_fd, 4096)
            except OSError:
                break
            if not data:
                break
            text = data.decode('utf-8', errors='replace')
            try:
                original_stderr.write(text)
                original_stderr.flush()
                _log_fh.write(text)
            except (OSError, ValueError):
                break
        try:
            os.close(r_fd)
        except OSError:
            pass

    _tee_thread = threading.Thread(target=_tee_fd2, daemon=True)
    _tee_thread.start()

    _term_fh = original_stderr
    sys.stdout = _Tee(sys.__stdout__, _log_fh)
    sys.stderr = _Tee(original_stderr, _log_fh)

    for level, text in _pending:
        if level >= _file_level:
            _log_fh.write(text + '\n')
    _pending.clear()

    step(f"Logging to {log_path}")
    return log_path


def log_dir() -> Path | None:
    return _logs_dir


def _shutdown() -> None:
    progress_end()
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.flush()
        except Exception:
            pass
    if _saved_stderr_fd is not None:
        try:
            os.dup2(_saved_stderr_fd, 2)
        except OSError:
            pass
    if _tee_thread is not None:
        _tee_thread.join(timeout=2.0)
    if _log_fh is not None:
        try:
            _log_fh.flush()
        except (OSError, ValueError):
            pass


atexit.register(_shutdown)


# ---------------------------------------------------------------------------
# Subprocesses
# ---------------------------------------------------------------------------

def tools_log_path() -> Path | None:
    if _logs_dir is None:
        return None
    return _logs_dir / 'tools.log'


def _tail(path: Path, n_lines: int = 20, max_bytes: int = 65536) -> list:
    try:
        size = path.stat().st_size
        with open(path, 'rb') as fh:
            if size > max_bytes:
                fh.seek(size - max_bytes)
            blob = fh.read()
    except OSError:
        return []
    text = blob.decode('utf-8', errors='replace').replace('\r', '\n')
    lines = [l.rstrip() for l in text.splitlines() if l.strip()]
    return lines[-n_lines:]


def run_tool(cmd, name: str, *, stdout=None, env=None, check: bool = True,
             cwd=None):
    argv = [str(c) for c in cmd]
    detail(f"  $ {shlex.join(argv)}")

    path = tools_log_path()
    if path is None:
        return subprocess.run(argv, check=check, stdout=stdout, env=env, cwd=cwd)

    with open(path, 'a', buffering=1) as fh:
        fh.write(
            f"\n# {'─' * 68}\n"
            f"# {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  {name}\n"
            f"# $ {shlex.join(argv)}\n"
        )
        fh.flush()
        try:
            return subprocess.run(
                argv, check=check,
                stdout=(fh if stdout is None else stdout), stderr=fh,
                env=env, cwd=cwd,
            )
        except subprocess.CalledProcessError as exc:
            fh.flush()
            error(f"{name} failed with exit status {exc.returncode}.")
            error(f"  command: {shlex.join(argv)}")
            error(f"  last lines of {path}:")
            for line in _tail(path):
                error(f"    {line[:200]}")
            raise

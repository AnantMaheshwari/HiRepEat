from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path


class ExternalToolError(RuntimeError):
    pass


TESORTER_DBS = (
    'rexdb',
    'rexdb-plant',
    'rexdb-metazoa',
    'rexdb-v3',
    'rexdb-plantv3',
    'rexdb-metazoav3',
    'rexdb-pnas',
    'rexdb-line',
    'gydb',
    'sine',
)

_IMAGE_OVERRIDE: Path | None = None

_NOT_FOUND = """TEsorter was not found, and classification is required."""


# ---------------------------------------------------------------------------
# Image / runtime discovery
# ---------------------------------------------------------------------------

def set_image(path) -> None:
    global _IMAGE_OVERRIDE
    _IMAGE_OVERRIDE = Path(path).resolve() if path else None


def image_path() -> Path | None:
    if _IMAGE_OVERRIDE is not None:
        return _IMAGE_OVERRIDE
    env = os.environ.get('HIREPEAT_IMAGE')
    return Path(env).resolve() if env else None


def in_container() -> bool:
    return bool(os.environ.get('HIREPEAT_IN_CONTAINER'))


def _container_runtime() -> str | None:
    return shutil.which('apptainer') or shutil.which('singularity')


def _search_path() -> str:
    """PATH to search for TEsorter, with $HIREPEAT_TOOL_PATH prepended."""
    base = os.environ.get('PATH', os.defpath)
    extra = os.environ.get('HIREPEAT_TOOL_PATH')
    return f"{extra}:{base}" if extra else base


def tool_env() -> dict | None:
    extra = os.environ.get('HIREPEAT_TOOL_PATH')
    if not extra:
        return None
    env = os.environ.copy()
    env['PATH'] = _search_path()
    return env


def _which(name: str) -> str | None:
    return shutil.which(name, path=_search_path())


def _which_tesorter() -> str | None:
    return _which('TEsorter')


def have_tesorter() -> bool:
    """True if TEsorter can be run at all, either way."""
    return bool(_which_tesorter()) or image_path() is not None


# ---------------------------------------------------------------------------
# Building the command
# ---------------------------------------------------------------------------

def tesorter_cmd(argv, bind_paths=()) -> list[str]:
    argv = [str(a) for a in argv]

    direct = _which_tesorter()
    if direct:
        return [direct] + argv

    sif = image_path()
    if sif is None:
        raise ExternalToolError(_NOT_FOUND)
    if not sif.exists():
        raise ExternalToolError(f"Container image does not exist: {sif}")

    runtime = _container_runtime()
    if runtime is None:
        raise ExternalToolError(
            f"Neither `apptainer` nor `singularity` is on PATH, so {sif} "
            "cannot be run."
        )

    cmd = [runtime, 'exec', '--cleanenv']
    for b in _bind_args(bind_paths):
        cmd += ['--bind', b]
    cmd.append(str(sif))
    return cmd + ['TEsorter'] + argv


def _bind_args(bind_paths) -> list[str]:
    dirs: list[Path] = []
    for p in bind_paths:
        if p is None:
            continue
        path = Path(p).resolve()
        d = path if path.is_dir() else path.parent
        if d not in dirs:
            dirs.append(d)

    # Drop any directory already covered by an ancestor in the list.
    minimal = []
    for d in dirs:
        if not any(d != other and other in d.parents for other in dirs):
            minimal.append(str(d))
    return minimal

def _run(cmd, env=None) -> str:
    """Run a probe command, returning stdout+stderr, or '' on any failure."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                           env=env)
    except (OSError, subprocess.SubprocessError):
        return ''
    return (r.stdout or '') + (r.stderr or '')


def _search(pattern, text, default=None):
    m = re.search(pattern, text)
    return m.group(1) if m else default


def probe_versions() -> dict:
    info: dict[str, object] = {
        'mode': None,
        'image': None,
        'tesorter_path': None,
        'tesorter': None,
        'hmmer': None,
        'blast': None,
    }

    sif = image_path()
    direct = _which_tesorter()
    env = None

    if direct:
        info['mode'] = 'container' if in_container() else 'path'
        info['tesorter_path'] = direct
        prefix: list[str] = []
        env = tool_env()
    elif sif is not None:
        info['mode'] = 'image'
        info['image'] = str(sif)
        runtime = _container_runtime()
        if runtime is None:
            return info
        prefix = [runtime, 'exec', '--cleanenv', str(sif)]
        labels = _sif_labels(sif)
        if labels:
            info['tesorter'] = labels.get('io.hirepeat.tesorter')
            info['hmmer'] = labels.get('io.hirepeat.hmmer')
            info['blast'] = labels.get('io.hirepeat.rmblast')
            info['hirepeat_image_version'] = labels.get(
                'org.opencontainers.image.version')
            info['hirepeat_image_revision'] = labels.get(
                'org.opencontainers.image.revision')
            if all(info[k] for k in ('tesorter', 'hmmer', 'blast')):
                return info
    else:
        return info

    def probe(tool, args, pattern):
        if prefix:
            cmd = prefix + [tool] + args
        else:
            exe = _which(tool)
            if exe is None:
                return None
            cmd = [exe] + args
        return _search(pattern, _run(cmd, env=env))

    if not info['tesorter']:
        info['tesorter'] = probe(
            'TEsorter', ['--version'], r'TEsorter\s+(\S+)')
    if not info['hmmer']:
        info['hmmer'] = probe(
            'hmmsearch', ['-h'], r'#\s*HMMER\s+(\S+)')
    if not info['blast']:
        info['blast'] = probe(
            'blastp', ['-version'], r'blastp:\s*(\S+)')

    return info


def _sif_labels(sif) -> dict:
    """Labels baked into a SIF, or {} if they cannot be read."""
    runtime = _container_runtime()
    if runtime is None:
        return {}
    out = _run([runtime, 'inspect', '--json', str(sif)])
    if not out.strip():
        return {}
    try:
        data = json.loads(out)
    except ValueError:
        return {}
    # Apptainer has moved this around between versions; accept both shapes.
    if isinstance(data.get('data'), dict):
        attrs = data['data'].get('attributes', {})
        labels = attrs.get('labels')
        if isinstance(labels, dict):
            return labels
    labels = data.get('attributes', {}).get('labels')
    return labels if isinstance(labels, dict) else {}


def preflight() -> dict:
    if not have_tesorter():
        raise ExternalToolError(_NOT_FOUND)
    return probe_versions()

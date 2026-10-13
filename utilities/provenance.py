"""Record what produced a run: parameters, inputs, and code version.
"""

import os
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from utilities import log

_REPO = Path(__file__).parent.parent

_INPUT_KEYS = ('genome', 'genome_list', 'genome_paths')


def _scalar(value):
    """Format one Python value as a YAML scalar."""
    if value is None:
        return 'null'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        return repr(value)
    text = str(value)
    return '"' + text.replace('\\', '\\\\').replace('"', '\\"') + '"'


def _git_state():
    def _git(*cmd):
        return subprocess.run(
            ('git',) + cmd, cwd=_REPO,
            capture_output=True, text=True, check=True,
        ).stdout
    try:
        commit = _git('rev-parse', '--short', 'HEAD').strip()
        dirty  = bool(_git('status', '--porcelain').strip())
        return commit, dirty
    except (OSError, subprocess.CalledProcessError):
        return None, None


def _version():
    """HiRepEat's own version.
    """
    env = os.environ.get('HIREPEAT_VERSION')
    if env:
        return env
    try:
        from importlib.metadata import version
        return version('HiRepEat')
    except Exception:
        return None


def render_run_params(args, output_dir, base, genome_bp=None,
                      tool_versions=None):
    """Render the run-parameters YAML as a string."""
    commit, dirty = _git_state()
    output_dir = Path(output_dir)

    lines = [
        '# repeat_framework run parameters',
        'run:',
        f'  started: {_scalar(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))}',
        f'  command: {_scalar(shlex.join([sys.executable] + sys.argv))}',
        f'  cwd: {_scalar(Path.cwd())}',
        f'  repo: {_scalar(_REPO)}',
        f'  version: {_scalar(_version())}',
        f'  git_commit: {_scalar(commit)}',
        f'  git_dirty: {_scalar(dirty)}',
        f'  python: {_scalar(sys.version.split()[0])}',
    ]

    if tool_versions:
        lines.append('tools:')
        for key in ('mode', 'image', 'tesorter_path', 'tesorter', 'hmmer',
                    'blast', 'hirepeat_image_version', 'hirepeat_image_revision'):
            if key in tool_versions and tool_versions[key] is not None:
                lines.append(f'  {key}: {_scalar(tool_versions[key])}')

    lines.append('input:')

    if args.genome:
        lines.append(f'  genome: {_scalar(Path(args.genome).resolve())}')
        if genome_bp is not None:
            lines.append(f'  genome_bp: {_scalar(genome_bp)}')
    if args.genome_list:
        lines.append(f'  genome_list: {_scalar(Path(args.genome_list).resolve())}')
        lines.append('  genomes:')
        for path in args.genome_paths:
            lines.append(f'    - {_scalar(Path(path).resolve())}')

    lines += [
        'output:',
        f'  dir: {_scalar(output_dir.resolve())}',
        f'  log: {_scalar((output_dir / f"{base}_pipeline.log").resolve())}',
        'params:',
    ]
    for key, value in sorted(vars(args).items()):
        if key in _INPUT_KEYS:
            continue
        lines.append(f'  {key}: {_scalar(value)}')

    lines.append('constants:')
    for key, value in sorted(_constants().items()):
        lines.append(f'  {key}: {_scalar(value)}')

    return '\n'.join(lines) + '\n'


def _constants():
    from engine import constants

    out = {
        name: getattr(constants, name)
        for name in vars(constants)
        if name.isupper() and not name.startswith('_')
    }

    out['TESORTER_GYDB'] = constants.tesorter_gydb()
    return out


def _unscalar(text):
    text = text.strip()
    if text == 'null':
        return None
    if text in ('true', 'false'):
        return text == 'true'
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        return text[1:-1].replace('\\"', '"').replace('\\\\', '\\')
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def read_run_params(path):

    path = Path(path)
    if not path.exists():
        return {}
    try:
        text = path.read_text()
    except OSError:
        return {}

    out: dict = {}
    section: dict | None = None
    last_key = None
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith('#'):
            continue
        indent = len(raw) - len(raw.lstrip(' '))
        line = raw.strip()

        if indent == 0:
            if line.endswith(':'):
                section = out.setdefault(line[:-1], {})
                last_key = None
            continue
        if section is None:
            continue

        if line.startswith('- '):
            if isinstance(section.get(last_key), list):
                section[last_key].append(_unscalar(line[2:]))
            continue

        key, sep, value = line.partition(':')
        if not sep:
            continue
        last_key = key.strip()
        section[last_key] = [] if not value.strip() else _unscalar(value)

    return out


def write_run_params(path, args, output_dir, base, genome_bp=None,
                     tool_versions=None):
    """Write run_params.yaml, and echo it into the log.  Returns the path."""
    text = render_run_params(args, output_dir, base, genome_bp=genome_bp,
                             tool_versions=tool_versions)
    path = Path(path)
    path.write_text(text)
    # The full block is a record, not news: it belongs in the log (and at -v),
    # while the terminal only needs to know where it landed.
    log.block(log.INFO, text)
    log.step(f"Run parameters written to {path}")
    return path

"""Checked links to public repository files, without copying private state."""
from pathlib import Path
import re
from urllib.parse import quote

from docutils import nodes, utils


def repo_file(name, rawtext, text, lineno, inliner, options=None, content=None):
    env = inliner.document.settings.env
    text = utils.unescape(text)
    match = re.fullmatch(r'(.+?)\s*<([^>]+)>', text)
    title, target = (match.group(1).strip(), match.group(2).strip()) if match else (text, text)
    path_text, _, fragment = target.partition('#')
    root = Path(env.app.config.repo_root).resolve()
    candidate = (root/path_text).resolve()
    forbidden = ('.local', '.git', '_build', '__pycache__')
    try:
        relative = candidate.relative_to(root)
        if (not candidate.is_file() or path_text.startswith('/') or
            any(part in forbidden for part in relative.parts) or
            relative.name.endswith('.local') or str(relative).startswith('artifacts/scalar-frontier/')):
            raise ValueError
    except ValueError:
        message = inliner.reporter.error(f'repo-file target is missing or private: {target}', line=lineno)
        return [inliner.problematic(rawtext, rawtext, message)], [message]
    env.note_dependency(str(candidate))
    url = f'{env.app.config.repo_url}/blob/{env.app.config.repo_revision}/{quote(relative.as_posix(), safe="/")}'
    if fragment:
        url += '#'+quote(fragment, safe='-')
    return [nodes.reference(rawtext, title, refuri=url, **(options or {}))], []


def setup(app):
    app.add_config_value('repo_root', '', 'env')
    app.add_config_value('repo_url', '', 'html')
    app.add_config_value('repo_revision', 'main', 'html')
    app.add_role('repo-file', repo_file)
    return {'version': '1.0', 'parallel_read_safe': True, 'parallel_write_safe': True}

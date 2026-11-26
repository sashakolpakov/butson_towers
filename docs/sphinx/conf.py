"""Sphinx build for the public construction and verification manual."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'src'), str(Path(__file__).parent/'_ext')]
project = 'Butson Towers'
author = 'Alexander Kolpakov'
copyright = '2026, Alexander Kolpakov'
version = release = '0.1.0'
needs_sphinx = '8.2'
extensions = [
    'myst_parser', 'sphinx.ext.autodoc', 'sphinx.ext.viewcode',
    'sphinx.ext.doctest', 'sphinx.ext.mathjax', 'project_links',
]
source_suffix = {'.md': 'markdown', '.rst': 'restructuredtext'}
root_doc = 'index'
exclude_patterns = ['_build', '_ext/**', '**/*.local', '**/.local/**', '.DS_Store']
myst_enable_extensions = ['dollarmath', 'amsmath', 'colon_fence', 'deflist']
myst_heading_anchors = 4
myst_dmath_double_inline = False
nitpicky = True
numfig = True
math_numfig = True
html_theme = 'furo'
html_title = 'Butson Towers: constructions, code and verification'
html_static_path = ['_static']
html_css_files = ['manual.css']
html_theme_options = {
    'source_repository': 'https://github.com/sashakolpakov/butson_towers/',
    'source_branch': 'main',
    'source_directory': 'docs/sphinx/',
}
html_show_sourcelink = True
html_last_updated_fmt = None
autodoc_member_order = 'bysource'
autodoc_typehints = 'none'
autodoc_preserve_defaults = True
viewcode_follow_imported_members = False
mathjax3_config = {
    'tex': {'tags': 'ams'},
}
doctest_global_setup = "import os; os.chdir(" + repr(str(ROOT)) + ")"
latex_engine = 'xelatex'
latex_documents = [('index', 'butson-towers-manual.tex',
                    'Butson Towers: Constructions, Code and Verification',
                    author, 'manual')]
latex_show_urls = 'no'
latex_elements = {
    'papersize': 'a4paper', 'pointsize': '10pt',
    'fontpkg': r'\setmainfont{Latin Modern Roman}\setsansfont{Latin Modern Sans}\setmonofont{Latin Modern Mono}',
}
repo_root = str(ROOT)
repo_url = 'https://github.com/sashakolpakov/butson_towers'
repo_revision = 'main'

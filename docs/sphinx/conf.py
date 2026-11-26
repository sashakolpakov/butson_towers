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
html_baseurl = 'https://sashakolpakov.github.io/butson_towers/'
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
add_module_names = False
python_maximum_signature_line_length = 88
viewcode_follow_imported_members = False
mathjax3_config = {
    'tex': {'tags': 'ams'},
}
doctest_global_setup = "import os; os.chdir(" + repr(str(ROOT)) + ")"
latex_engine = 'xelatex'
latex_use_xindy = False
latex_documents = [('index', 'butson-towers-manual.tex',
                    'Butson Towers: Constructions, Code and Verification',
                    author, 'manual')]
latex_show_urls = 'no'
latex_elements = {
    'papersize': 'a4paper', 'pointsize': '10pt',
    'preamble': r'''
\makeatletter
\renewenvironment{sphinxtheindex}{%
  \clearpage\chapter*{\indexname}\phantomsection
  \addcontentsline{toc}{chapter}{\indexname}%
  \markboth{\MakeUppercase\indexname}{\MakeUppercase\indexname}%
  \thispagestyle{plain}\small\raggedright
  \parindent\z@\parskip\z@\@plus .3\p@\relax
  \let\item\@idxitem
}{\clearpage}
\makeatother
''',
    'fontpkg': r'''
\setmainfont{lmroman10-regular.otf}[BoldFont=lmroman10-bold.otf,ItalicFont=lmroman10-italic.otf,BoldItalicFont=lmroman10-bolditalic.otf,SmallCapsFont=lmromancaps10-regular.otf,FontFace={m}{scit}{Font=lmromancaps10-oblique.otf}]
\setsansfont{lmsans10-regular.otf}[BoldFont=lmsans10-bold.otf,ItalicFont=lmsans10-oblique.otf,BoldItalicFont=lmsans10-boldoblique.otf]
\setmonofont{DejaVuSansMono.ttf}[BoldFont=DejaVuSansMono-Bold.ttf,ItalicFont=DejaVuSansMono-Oblique.ttf,BoldItalicFont=DejaVuSansMono-BoldOblique.ttf,Scale=MatchLowercase]
''',
}
repo_root = str(ROOT)
repo_url = 'https://github.com/sashakolpakov/butson_towers'
repo_revision = 'main'

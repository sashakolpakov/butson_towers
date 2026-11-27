#!/usr/bin/env python3
"""Build and bind version2 of existing Hexagon article2611.01070, offline."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TITLE = 'Butson towers and high-rank Hadamard matrices'
IDENTIFIER = '2611.01070'
LEAN_COMMIT = '6e515aecc448181c52d8c1eba3f3b0db0d9d5f11'
LEAN_PATH = 'release/additional_source_verification/butson-towers-additional-sources.zip'
ADMIN_FIELDS = {'authors', 'submittingContributorIndex', 'sponsorName',
                'submittingOrganizationId', 'license', 'authorshipMode',
                'humanUnderstanding', 'provenance', 'correctnessConfirmed',
                'communityInterestConfirmed', 'citationsConfirmed', 'rightsConfirmed'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def prose(value):
    value = re.sub(r'(?<!\\)%[^\n]*', '', value).replace('$', '').replace('~', ' ')
    require('\\' not in value and '{' not in value and '}' not in value,
            'Extend the metadata TeX parser for this changed title/abstract')
    return re.sub(r'\s+', ' ', value).strip()


def private_path(name):
    return any(p.lower() == '.local' or p.lower().endswith('.local') for p in Path(name).parts)


def source_bytes(relative):
    require(not private_path(relative), 'Private source path refused')
    path = ROOT/relative
    require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(ROOT),
            f'Missing or unsafe source: {relative}')
    require(all(not p.is_symlink() for p in path.parents if p.is_relative_to(ROOT)),
            f'Symlinked source directory: {relative}')
    return path.read_bytes()


def lean_pointer():
    package = ROOT/'release/additional_source_verification'
    report_bytes = (package/'verification.json').read_bytes()
    report = json.loads(report_bytes)
    archive_bytes = (ROOT/LEAN_PATH).read_bytes()
    require(report['success'] and report['proof_theorems_checked'] == 143
            and not report['project_axioms_or_placeholders_used'], 'Incomplete Lean audit')
    require(report['archive_sha256'] == digest(archive_bytes), 'Lean archive digest mismatch')
    manifest = json.loads((package/'source/MANIFEST.json').read_text())
    matched = 0
    # The formal tree and literal inputs must still equal the pinned source packet.
    for entry in manifest['files']:
        relative = entry['path']
        if relative.startswith('formal/') or relative.startswith('experiments/') or relative.startswith('scripts/') or relative.startswith('artifacts/'):
            canonical = ROOT/relative
            if canonical.is_file():
                require(digest(canonical.read_bytes()) == entry['sha256'],
                        f'Current formal dependency differs from the pinned packet: {relative}')
                matched += 1
    require(matched >= 70, 'The current formal closure was not compared')
    (HERE/'lean-verification.json').write_bytes(report_bytes)
    pointer = {'schema': 'butson-current-Lean-package-pointer-v1',
               'scope': '34 announced existence orders and nine all-parameter families;143 selected axiom reports',
               'proof_commit': LEAN_COMMIT,
               'source_package_url': f'https://github.com/sashakolpakov/butson_towers/raw/{LEAN_COMMIT}/{LEAN_PATH}',
               'source_package_sha256': digest(archive_bytes),
               'verification_record_sha256': digest(report_bytes),
               'payload_manifest_sha256': report['payload_manifest_sha256'],
               'current_formal_and_literal_source_files_matched': matched,
               'lean_toolchain': report['lean_toolchain'], 'mathlib_revision': report['mathlib_revision'],
               'proof_theorems_checked': 143,
               'finite_field_rank_comparison_in_Lean': False,
               'replay': 'Extract the source-only ZIP and run python3 verify.py; see itsREADME for pinned dependency/cache options.'}
    write_json(HERE/'proof-sources.json', pointer)
    return pointer


def metadata(source, pages, private_root):
    text = source.decode()
    title = prose(re.search(r'\\title\{([^{}]*)\}', text).group(1))
    abstract = prose(text.split(r'\begin{abstract}', 1)[1].split(r'\end{abstract}', 1)[0])
    require(title == TITLE and 'thirty-four new orders' in abstract and 'nine' in abstract,
            'Title or abstract count differs from the current theorem')
    parameter_path = private_root/'hexagon_parameters.local'
    require(parameter_path.stat().st_mode & 0o777 == 0o600, 'Private preferences must have mode0600')
    preferences = json.loads(parameter_path.read_text())
    require(preferences['schema'] == 'hexagon-fixed-submission-parameters-v1', 'Unknown preference format')
    defaults = preferences['fixed_administrative_parameters']
    require(set(defaults) == ADMIN_FIELDS, 'Unexpected fixed administrative fields')
    result = dict(defaults)
    result['authors'] = [{k: v for k, v in author.items() if k not in ('contributorKind', 'organizationId')}
                         for author in defaults['authors']]
    result.update(title=title, abstract=abstract, subject='math.CO', crossListSubjects=[],
                  mscCodes=['05B20', '15B34'],
                  tags=['Butson matrices', 'Hadamard matrices', 'computer-assisted proof',
                        'finite fields', 'restriction of scalars'],
                  comments=(f'{pages} pages. Revision of2611.01070v1:34 new orders and nine new infinite families '
                            '(previously18 and six). Adds a separate table of31 recovered known orders and an '
                            'exact proof that the order885924 sixth-root representative is inequivalent to '
                            'the standard Paley matrix. All34 announced orders and nine families have complete '
                            'Lean existence proofs from literal finite witnesses;143 selected axiom reports '
                            'and the pinned source-only proof package are provided. The finite-field rank '
                            'comparison has a separate exact executable certificate. Includes canonical '
                            'Python source recipes, literal data and standalone archive replays.'))
    result['provenance'] = defaults['provenance'].strip()+'\n\n'+(
        'The search implementation, exact arithmetic checks, automated cross-checks and initial text '
        'were developed with OpenAI Codex and automated collaborators. The paper distinguishes '
        'Lean existence proofs, executable finite certificates and literature novelty labels.')
    # Version comments are plain prose, including spaces around every numeric reference.
    result['comments'] = (result['comments'].replace('of2611', 'of 2611').replace('v1:34', 'v1: 34')
                          .replace('previously18', 'previously 18').replace('of31', 'of 31')
                          .replace('order885924', 'order-885924').replace('All34', 'All 34')
                          .replace(';143', '; 143'))
    write_json(HERE/'metadata.template.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--private-root', type=Path, default=ROOT)
    args = parser.parse_args()
    if (HERE/'submission-status.json').exists():
        manifest = json.loads((HERE/'submission-manifest.json').read_text())
        require((HERE/'main.tex').read_bytes() == source_bytes('manuscript/hadamard_7572_announcement.tex'),
                'This submitted package is frozen; prepare a later manuscript version in a new release directory')
        for name, field in [('source.zip', 'source_archive_sha256'),
                            ('announcement.pdf', 'reference_pdf_sha256'),
                            ('metadata.template.json', 'metadata_sha256')]:
            require(digest((HERE/name).read_bytes()) == manifest[field], 'A frozen release digest differs')
        print(json.dumps({'submitted_release_frozen': True, 'identifier': IDENTIFIER,
                          'version': IDENTIFIER+'v2', 'all_bound_release_hashes_verified': True,
                          'files_rewritten': False}))
        return
    payload = json.loads((HERE/'payload.json').read_text())
    require(payload['schema'] == 'butson-hexagon-v2-payload-v1', 'Unexpected payload manifest')
    source = source_bytes('manuscript/hadamard_7572_announcement.tex')
    (HERE/'main.tex').write_bytes(source)
    command = ['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', '-outdir=build', 'main.tex']
    built = subprocess.run(command, cwd=HERE, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (HERE/'build/latexmk-output.log').write_text(built.stdout)
    require(built.returncode == 0, 'LaTeX failed; inspect the retained local build log')
    log = (HERE/'build/main.log').read_text()
    warnings = [line for line in log.splitlines() if 'Warning:' in line or 'Overfull' in line or 'Underfull' in line]
    require(not warnings, 'Resolve all manuscript typesetting warnings before packaging')
    page_match = re.search(r'Output written on .*?\((\d+) pages?', log, re.S)
    require(page_match is not None, 'PDF page count missing')
    pages = int(page_match[1])
    pdf = (HERE/'build/main.pdf').read_bytes()
    (HERE/'announcement.pdf').write_bytes(pdf)
    write_json(HERE/'build.json', {'source': 'main.tex', 'source_sha256': digest(source),
                'reference_pdf': 'announcement.pdf', 'pdf_sha256': digest(pdf), 'pages': pages,
                'command': command, 'warnings': warnings,
                'source_equals_current_manuscript': True, 'all_pages_visually_reviewed': False})
    meta = metadata(source, pages, args.private_root.resolve())
    pointer = lean_pointer()
    (HERE/'manuscript_data.py').write_bytes(source_bytes('scripts/verify_manuscript_data.py'))
    files = {'main.tex': source}
    for relative in payload['repository_files']:
        require(Path(relative).suffix in ('.py', '.json'), f'Unsupported ancillary extension: {relative}')
        files['anc/'+relative] = source_bytes(relative)
    for relative, destination in payload['release_files'].items():
        require(not private_path(relative) and destination.startswith('anc/'), 'Unsafe release member')
        require(Path(relative).suffix in ('.py', '.json'), 'Unsupported release file extension')
        files[destination] = (HERE/relative).read_bytes()
    require(len(files) == 1+len(payload['repository_files'])+len(payload['release_files']), 'Archive path collision')
    ancillary_manifest = {'schema': 'butson-hexagon-v2-ancillary-v1',
                          'purpose': 'Canonical source recipes and exact verification for34 new orders,nine families,and31 known recovered orders',
                          'included_files': {name: digest(value) for name, value in sorted(files.items())},
                          'new_orders': 34, 'new_families': 9, 'known_orders': 31,
                          'Lean_axiom_reports': 143, 'rank_comparison_Lean_checked': False,
                          'dependencies': payload['dependencies'],
                          'replays': payload['replays']}
    files['anc/manifest.json'] = (json.dumps(ancillary_manifest, indent=2)+'\n').encode()
    secret = (args.private_root/'hexagon_key.local').read_bytes().strip()
    preference_bytes = (args.private_root/'hexagon_parameters.local').read_bytes()
    require(secret, 'Credential file is empty')
    for name, value in files.items():
        require(not private_path(name) and not Path(name).is_absolute() and '..' not in Path(name).parts,
                'Unsafe archive member name')
        require(len(value) <= 4194304, 'Text source exceeds the service per-file limit')
        require(secret not in value and preference_bytes not in value, 'Private data found in archive')
    archive_path = HERE/'source.zip'
    with zipfile.ZipFile(archive_path, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(files.items()):
            entry = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
            entry.external_attr = 0o100644 << 16
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, data)
    manifest = {'schema': 'butson-hexagon-existing-article-version-v2', 'identifier': IDENTIFIER,
                'previous_version': IDENTIFIER+'v1', 'expected_version': IDENTIFIER+'v2',
                'source_archive': 'source.zip', 'source_archive_sha256': digest(archive_path.read_bytes()),
                'reference_pdf_sha256': digest(pdf), 'metadata_sha256': digest((HERE/'metadata.template.json').read_bytes()),
                'main_source_sha256': digest(source), 'pages': pages,
                'new_orders': 34, 'new_families': 9, 'known_orders_recovered': 31,
                'distinct_finite_completed_orders': 65, 'Lean_selected_axiom_reports': 143,
                'pinned_Lean_package_sha256': pointer['source_package_sha256'],
                'source_archive_members': [{'path': name, 'bytes': len(value), 'sha256': digest(value)}
                                           for name, value in sorted(files.items())],
                'private_local_files_excluded': True, 'network_requests_performed_by_builder': False}
    write_json(HERE/'submission-manifest.json', manifest)
    write_json(HERE/'version-target.json', {'identifier': IDENTIFIER, 'previous_version': IDENTIFIER+'v1',
               'expected_version': IDENTIFIER+'v2', 'title': TITLE,
               'new_article_creation_permitted': False,
               'operation': 'POST /api/v1/submissions/{originalDraftId}/revise'})
    print(json.dumps({'built': True, 'pages': pages, 'archive_members': len(files),
                      'archive_bytes': archive_path.stat().st_size,
                      'abstract_words': len(meta['abstract'].split()), 'network_requests': False}))


if __name__ == '__main__':
    main()

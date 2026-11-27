#!/usr/bin/env python3
"""Upload only a revision of existing Hexagon article 2611.01070.

No request in this client can create a new article. Credentials and all
private draft receipts remain under ignored local paths. A successful
preview and exact expanded-source comparison are required before commit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from datetime import datetime, timezone
from io import BytesIO
from urllib.error import HTTPError
from urllib.request import HTTPRedirectHandler, Request, build_opener

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ORIGIN = 'https://hexagonmath.org'
IDENTIFIER = '2611.01070'
PREVIOUS_VERSION = IDENTIFIER + 'v1'
EXPECTED_VERSION = IDENTIFIER + 'v2'
TITLE = 'Butson towers and high-rank Hadamard matrices'
CREDENTIAL_ROOT = ROOT
PRIVATE = ROOT / '.local' / 'hexagon-2611-01070-v2'
PREVIOUS_STATE = ROOT / '.local' / 'hexagon-butson-towers-v3' / 'state.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def private_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(path.parent, 0o700)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(descriptor, 'w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
    os.chmod(path, 0o600)


def credential():
    path = CREDENTIAL_ROOT / 'hexagon_key.local'
    require(path.stat().st_mode & 0o777 == 0o600, 'Credential file must have mode0600')
    raw = path.read_text().strip()
    match = re.fullmatch(r'(?:export\s+)?(?:HEXAGON_API_KEY|HEXAGON_TOKEN)\s*=\s*(.+)', raw, re.S)
    value = match.group(1).strip() if match else raw
    if value[:1] == value[-1:] and value[:1] in ('\"', "'"):
        value = value[1:-1].strip()
    if value.lower().startswith('bearer '):
        value = value[7:].strip()
    require(value and '\n' not in value and '\r' not in value, 'Invalid credential format')
    return value


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError('Authenticated redirects are refused')


class Session:
    def __init__(self):
        self.token = credential()
        self.opener = build_opener(NoRedirect())
        PRIVATE.mkdir(parents=True, exist_ok=True, mode=0o700)
        os.chmod(PRIVATE, 0o700)

    def request(self, method, path, body=None, raw=None, binary=False):
        require(path.startswith('/api/v1/'), 'Only the official API is allowed')
        require(not (method == 'POST' and path.rstrip('/') == '/api/v1/submissions'),
                'Creating a new article is prohibited')
        require('://' not in path and not any(c in path for c in ('\n', '\r')), 'Invalid API path')
        require(body is None or raw is None, 'Ambiguous request body')
        data = json.dumps(body).encode() if body is not None else raw
        if method == 'POST' and data is None:
            data = b''
        headers = {'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'}
        if path.startswith('/api/v1/submissions/'):
            headers['Authorization'] = 'Bearer ' + self.token
        else:
            require(method == 'GET' and path.startswith('/api/v1/works/'),
                    'This client allows only public work reads and the submission revision API')
        if data is not None:
            headers['Content-Type'] = 'application/octet-stream' if raw is not None else 'application/json'
        request = Request(ORIGIN+path, data=data, headers=headers, method=method)
        try:
            with self.opener.open(request, timeout=60) as response:
                content = response.read()
                value = content if binary else (json.loads(content) if content else {})
                status = response.status
        except HTTPError as error:
            content = error.read()
            private_json(PRIVATE/'last-error.json', {'method': method, 'path': path,
                         'status': error.code, 'body': content.decode(errors='replace')})
            raise ValueError(f'Hexagon returned HTTP{error.code}; response saved privately') from None
        if binary:
            return value
        records = sorted(PRIVATE.glob('response-*.json'))
        private_json(PRIVATE/f'response-{len(records)+1:03d}.json',
                     {'method': method, 'path': path, 'status': status, 'body': value})
        return value


def previous_target():
    previous = json.loads(PREVIOUS_STATE.read_text())
    require(previous['identifier'] == IDENTIFIER and previous['version_id'] == PREVIOUS_VERSION
            and previous['committed'] is True, 'Prior receipt does not match the requested existing article')
    require(re.fullmatch(r'[A-Za-z0-9_-]{1,128}', previous['draft_id']), 'Invalid prior draft id')
    return previous['draft_id']


def inspect(session):
    old = session.request('GET', '/api/v1/submissions/'+previous_target())
    public = session.request('GET', '/api/v1/works/'+PREVIOUS_VERSION)['work']
    require(public['identifier'] == IDENTIFIER and public['versionId'] == PREVIOUS_VERSION
            and public['title'] == TITLE and public['status'] == 'published', 'Public prior article mismatch')
    private_json(PRIVATE/'prior-public-work.json', public)
    print(json.dumps({'existing_article_verified': True, 'identifier': IDENTIFIER,
                      'previous_version': PREVIOUS_VERSION, 'public_status': public['status'],
                      'prior_deposit_status': old.get('status'), 'private_deposit_fields': list(old)}))



def local_release():
    manifest = json.loads((HERE/'submission-manifest.json').read_text())
    metadata_bytes = (HERE/'metadata.template.json').read_bytes()
    metadata = json.loads(metadata_bytes)
    require(manifest['identifier'] == IDENTIFIER and manifest['previous_version'] == PREVIOUS_VERSION
            and manifest['expected_version'] == EXPECTED_VERSION and metadata['title'] == TITLE,
            'The local release does not target the requested article version')
    require(sha(metadata_bytes) == manifest['metadata_sha256']
            and sha((HERE/'source.zip').read_bytes()) == manifest['source_archive_sha256']
            and sha((HERE/'announcement.pdf').read_bytes()) == manifest['reference_pdf_sha256'],
            'Local release hashes changed')
    replay = json.loads((HERE/'archive-verification.json').read_text())
    require(replay['success'] and replay['source_archive_sha256'] == manifest['source_archive_sha256']
            and not replay['check_only'] and len(replay['replays']) >= 12
            and all(item['success'] for item in replay['replays']),
            'Fresh extracted source replays must pass before upload')
    return manifest, metadata


def state_for(manifest):
    path = PRIVATE/'state.json'
    if not path.exists():
        return None
    state = json.loads(path.read_text())
    require(state['identifier'] == IDENTIFIER and state['expected_version'] == EXPECTED_VERSION
            and state['source_archive_sha256'] == manifest['source_archive_sha256']
            and state['metadata_sha256'] == manifest['metadata_sha256'],
            'The retained revision state belongs to different release bytes')
    return state


def save_state(state):
    private_json(PRIVATE/'state.json', state)


def verify_draft_target(value):
    require(value.get('workIdentifier') == IDENTIFIER,
            'Revision draft does not belong to the existing requested article')
    version = value.get('versionId')
    # Hexagon retains the base version on an editable revision until commit.
    # The distinct draft was created only through that base deposit's revise route;
    # the final receipt must identify the exact new version.
    base_pointer = (version == PREVIOUS_VERSION and value.get('status') == 'uploading'
                    and not value.get('receipt'))
    require(version in (None, EXPECTED_VERSION) or base_pointer,
            'Revision draft targets an unexpected version')


def verify_remote_files(files, manifest):
    expected = {item['path']: item for item in manifest['source_archive_members']}
    require(len(files) == len(expected) and {item['path'] for item in files} == set(expected),
            'The expanded remote source tree differs from the local archive')
    for item in files:
        local = expected[item['path']]
        require(item['sha256'] == local['sha256'] and item['sizeBytes'] == local['bytes'],
                'A remote expanded file has different bytes')



def clear_inherited_tree(session):
    manifest, _ = local_release()
    state = state_for(manifest)
    require(state is not None and state.get('created_by_revision_of_original_deposit') is True
            and state['draft_id'] != previous_target() and not state.get('committed')
            and 'upload_id' not in state,
            'Only inherited files of the new editable revision may be cleared')
    base = '/api/v1/submissions/'+state['draft_id']
    draft = session.request('GET', base)
    verify_draft_target(draft)
    require(draft['status'] == 'uploading', 'The revision is not editable')
    parent = session.request('GET', '/api/v1/submissions/'+previous_target())
    require(parent['status'] == 'published' and parent['workIdentifier'] == IDENTIFIER
            and parent['versionId'] == PREVIOUS_VERSION, 'The published parent target changed')
    prior = {item['relativePath']: (item['sha256'], item['sizeBytes'], item['kind'])
             for item in parent['uploads']}
    current = draft['uploads']
    require(all(item['relativePath'] in prior
                and (item['sha256'], item['sizeBytes'], item['kind']) == prior[item['relativePath']]
                for item in current), 'An editable file is not an unchanged clone of version1')
    private_json(PRIVATE/'inherited-source-tree.json', current)
    # Expanded members belong to their archive and are removed with it.
    for item in current:
        if item['kind'] == 'source-archive':
            session.request('DELETE', base+'/uploads/'+item['id'])
    remaining = session.request('GET', base)['uploads']
    for item in remaining:
        require(item['relativePath'] in prior
                and (item['sha256'], item['sizeBytes'], item['kind']) == prior[item['relativePath']],
                'An editable file changed while clearing inherited copies')
        session.request('DELETE', base+'/uploads/'+item['id'])
    count = len(current)
    refreshed = session.request('GET', base)
    require(not refreshed['uploads'], 'Inherited file copies remain in the revision draft')
    state['inherited_v1_file_copies_removed'] = True
    state['inherited_v1_file_copy_count'] = len(prior)
    save_state(state)
    print(json.dumps({'inherited_revision_copies_removed': count,
                      'published_parent_unchanged': True, 'identifier': IDENTIFIER}))

def stage(session):
    manifest, metadata = local_release()
    state = state_for(manifest)
    if state is None:
        inspect(session)
        original = previous_target()
        revised = session.request('POST', '/api/v1/submissions/'+original+'/revise')
        new_id = revised.get('draftId') or revised.get('id')
        if not new_id and isinstance(revised.get('draft'), dict):
            new_id = revised['draft'].get('id')
        require(isinstance(new_id, str) and re.fullmatch(r'[A-Za-z0-9_-]{1,128}', new_id)
                and new_id != original, 'The revision operation did not return a distinct editable draft')
        state = {'identifier': IDENTIFIER, 'previous_version': PREVIOUS_VERSION,
                 'expected_version': EXPECTED_VERSION, 'draft_id': new_id,
                 'created_by_revision_of_original_deposit': True,
                 'source_archive_sha256': manifest['source_archive_sha256'],
                 'metadata_sha256': manifest['metadata_sha256'], 'committed': False}
        save_state(state)
    require(not state['committed'], 'This revision has already been committed; read its status')
    base = '/api/v1/submissions/'+state['draft_id']
    draft = session.request('GET', base)
    verify_draft_target(draft)
    session.request('PATCH', base, body=metadata)
    if 'upload_id' not in state:
        upload = session.request('POST', base+'/uploads', body={
            'kind': 'source-archive', 'filename': 'butson-towers-2611.01070-v2.zip',
            'relativePath': 'butson-towers-2611.01070-v2.zip',
            'contentType': 'application/zip', 'sizeBytes': (HERE/'source.zip').stat().st_size})
        state.update(upload_id=upload['uploadId'], part_size=upload['partSizeBytes'])
        save_state(state)
    upload_base = base+'/uploads/'+state['upload_id']
    data = (HERE/'source.zip').read_bytes()
    if not state.get('upload_complete'):
        for index, offset in enumerate(range(0, len(data), state['part_size']), 1):
            session.request('PUT', upload_base+'/parts/'+str(index),
                            raw=data[offset:offset+state['part_size']])
        session.request('POST', upload_base+'/complete')
        state['upload_complete'] = True
        save_state(state)
    expanded = session.request('POST', upload_base+'/expand')
    verify_remote_files(expanded['files'], manifest)
    state['all_remote_expanded_file_hashes_match'] = True
    state['remote_expanded_file_count'] = len(expanded['files'])
    save_state(state)
    session.request('POST', upload_base+'/entry-point', body={'relativePath': 'main.tex'})
    static = session.request('GET', base+'/preview/check')
    require(static.get('ready') is True, 'Remote source pre-check is not ready')
    preview = session.request('POST', base+'/preview')
    state['preview_requested'] = True
    state['preview_digest'] = preview['digest']
    save_state(state)
    print(json.dumps({'revision_staged': True, 'identifier': IDENTIFIER,
                      'expected_version': EXPECTED_VERSION,
                      'remote_expanded_files_verified': len(expanded['files']),
                      'preview_status': (preview.get('preview') or {}).get('status')}))


def preview_status(session):
    manifest, _ = local_release()
    state = state_for(manifest)
    require(state is not None, 'No revision draft is staged')
    base = '/api/v1/submissions/'+state['draft_id']
    value = session.request('GET', base+'/preview')
    preview = value.get('preview') or {}
    print(json.dumps({'identifier': IDENTIFIER, 'expected_version': EXPECTED_VERSION,
                      'preview_status': preview.get('status'), 'current': preview.get('current'),
                      'gate': (value.get('gate') or {}).get('state')}))
    return value


def validate_receipt(value):
    receipt = value.get('receipt') or {}
    require(value.get('identifier') in (IDENTIFIER, 'hexagon:'+IDENTIFIER)
            and value.get('versionId') == EXPECTED_VERSION
            and receipt.get('identifier') == IDENTIFIER and receipt.get('versionId') == EXPECTED_VERSION
            and receipt.get('title') == TITLE,
            'The service receipt does not identify exactly version2 of the existing article')
    return receipt


def commit(session):
    from pypdf import PdfReader
    manifest, metadata = local_release()
    state = state_for(manifest)
    require(state is not None, 'No revision draft is staged')
    if state.get('committed'):
        print(json.dumps({'already_committed': True, 'identifier': IDENTIFIER,
                          'version': EXPECTED_VERSION}))
        return
    require(state.get('created_by_revision_of_original_deposit') is True
            and state.get('all_remote_expanded_file_hashes_match') is True,
            'The existing-article revision path and expanded source checks are required')
    base = '/api/v1/submissions/'+state['draft_id']
    draft = session.request('GET', base)
    verify_draft_target(draft)
    remote_metadata = draft['metadataJson']
    if isinstance(remote_metadata, str):
        remote_metadata = json.loads(remote_metadata)
    for name in ('title', 'abstract', 'comments', 'subject', 'mscCodes', 'license',
                 'authorshipMode', 'humanUnderstanding', 'provenance'):
        require(remote_metadata.get(name) == metadata.get(name), 'Remote revision metadata differs')
    preview = session.request('GET', base+'/preview')
    result = preview.get('preview') or {}
    require(result.get('status') == 'succeeded' and result.get('current') is True
            and result.get('hasPdf') is True and (preview.get('gate') or {}).get('allowed') is True
            and preview.get('digest') == state.get('preview_digest') == result.get('digest'),
            'A successful current server preview is required before commit')
    pdf = session.request('GET', base+'/preview/pdf', binary=True)
    parsed = PdfReader(BytesIO(pdf))
    require(len(parsed.pages) == manifest['pages'], 'Server preview page count differs')
    texts = [page.extract_text() for page in parsed.pages]
    require('thirty-four' in texts[0] and 'nine' in texts[0]
            and any('246090' in text and '442962' in text for text in texts),
            'Server preview lacks the revised abstract or rank comparison')
    descriptor = os.open(PRIVATE/'preview.pdf', os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(descriptor, 'wb') as stream:
        stream.write(pdf)
    state.update(server_preview_succeeded=True, server_preview_pages=len(parsed.pages),
                 server_preview_sha256=sha(pdf))
    save_state(state)
    committed = session.request('POST', base+'/commit')
    receipt = validate_receipt(committed)
    state.update(committed=True, version_id=EXPECTED_VERSION, receipt=receipt)
    save_state(state)
    private_json(PRIVATE/'receipt.local', committed)
    status(session)



def recover_receipt(session):
    manifest, _ = local_release()
    state = state_for(manifest)
    require(state is not None and state.get('server_preview_succeeded') is True,
            'No verified preview precedes the retained commit response')
    path = '/api/v1/submissions/'+state['draft_id']+'/commit'
    candidates = [json.loads(file.read_text()) for file in sorted(PRIVATE.glob('response-*.json'))]
    matches = [value for value in candidates if value['method'] == 'POST'
               and value['path'] == path and value['status'] == 202]
    require(matches, 'No successful commit response was retained')
    body = matches[-1]['body']
    receipt = validate_receipt(body)
    current = session.request('GET', '/api/v1/submissions/'+state['draft_id'])
    verify_draft_target(current)
    require(current.get('receipt') == receipt, 'Current service receipt differs from the accepted commit')
    state.update(committed=True, version_id=EXPECTED_VERSION, receipt=receipt)
    save_state(state)
    private_json(PRIVATE/'receipt.local', body)
    status(session)

def status(session):
    manifest, _ = local_release()
    state = state_for(manifest)
    require(state is not None, 'No revision exists')
    value = session.request('GET', '/api/v1/submissions/'+state['draft_id'])
    verify_draft_target(value)
    if state.get('committed'):
        receipt = state.get('receipt') or value.get('receipt')
        require(receipt['identifier'] == IDENTIFIER and receipt['versionId'] == EXPECTED_VERSION,
                'Retained receipt target mismatch')
        public = {'schema': 'hexagon-existing-article-version-receipt-v1',
                  'service': ORIGIN, 'identifier': IDENTIFIER, 'previous_version': PREVIOUS_VERSION,
                  'version_id': EXPECTED_VERSION, 'version_number': 2, 'title': TITLE,
                  'operation': 'revision_of_existing_published_article',
                  'new_article_created': False, 'status': value['status'],
                  'received_at': receipt['receivedAt'],
                  'checked_at_utc': datetime.now(timezone.utc).isoformat(),
                  'source_archive_sha256': manifest['source_archive_sha256'],
                  'metadata_sha256': manifest['metadata_sha256'],
                  'source_tree_sha256': receipt['sourceTreeSha256'],
                  'expanded_source_files': state['remote_expanded_file_count'],
                  'all_remote_expanded_file_hashes_match': True,
                  'server_preview_compile': 'succeeded', 'server_preview_pages': state['server_preview_pages'],
                  'public_record_url': ORIGIN+'/'+EXPECTED_VERSION if value['status'] == 'published' else None,
                  'processing_jobs': [{'kind': job.get('kind'), 'status': job.get('status')}
                                      for job in value.get('jobs', [])],
                  'private_receipt_saved_in_ignored_local_directory': True}
        (HERE/'submission-status.json').write_text(json.dumps(public, indent=2)+'\n')
    print(json.dumps({'identifier': IDENTIFIER, 'version': EXPECTED_VERSION,
                      'status': value['status'], 'revision_committed': bool(state.get('committed')),
                      'jobs': [{'kind': job.get('kind'), 'status': job.get('status')}
                               for job in value.get('jobs', [])]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=['inspect', 'clear-inherited', 'stage', 'preview', 'commit', 'recover-receipt', 'status'])
    parser.add_argument('--private-root', type=Path, default=ROOT,
                        help='checkout containing the ignored credential and prior receipt')
    args = parser.parse_args()
    global CREDENTIAL_ROOT, PRIVATE, PREVIOUS_STATE
    CREDENTIAL_ROOT = args.private_root.resolve()
    PRIVATE = CREDENTIAL_ROOT / '.local' / 'hexagon-2611-01070-v2'
    PREVIOUS_STATE = CREDENTIAL_ROOT / '.local' / 'hexagon-butson-towers-v3' / 'state.json'
    operation = {'inspect': inspect, 'clear-inherited': clear_inherited_tree,
                 'stage': stage, 'preview': preview_status,
                 'commit': commit, 'recover-receipt': recover_receipt, 'status': status}[args.phase]
    operation(Session())


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(1)

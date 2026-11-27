"""Prevent the revision uploader from creating a new article or changing target."""
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1]/'release/hexagon_2611_01070_v2/upload_version.py'
SPEC = importlib.util.spec_from_file_location('hexagon_existing_version', SCRIPT)
uploader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(uploader)


class FakeResponse:
    status = 200

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self):
        return b'{}'


class CapturingOpener:
    def __init__(self):
        self.requests = []

    def open(self, request, timeout):
        self.requests.append(request)
        return FakeResponse()


def session(monkeypatch, tmp_path):
    monkeypatch.setattr(uploader, 'PRIVATE', tmp_path)
    value = uploader.Session.__new__(uploader.Session)
    value.token = 'test-credential'
    value.opener = CapturingOpener()
    return value


def test_new_article_creation_is_refused_before_network(monkeypatch, tmp_path):
    value = session(monkeypatch, tmp_path)
    for path in ('/api/v1/submissions', '/api/v1/submissions/'):
        with pytest.raises(ValueError, match='Creating a new article is prohibited'):
            value.request('POST', path, body={'title': uploader.TITLE})
    assert not value.opener.requests


def test_authentication_stays_on_submission_api(monkeypatch, tmp_path):
    value = session(monkeypatch, tmp_path)
    value.request('GET', '/api/v1/works/2611.01070v1')
    assert not value.opener.requests[-1].has_header('Authorization')
    value.request('GET', '/api/v1/submissions/existing-draft')
    assert value.opener.requests[-1].get_header('Authorization') == 'Bearer test-credential'
    value.request('POST', '/api/v1/submissions/existing-draft/revise')
    assert value.opener.requests[-1].get_method() == 'POST'
    with pytest.raises(ValueError):
        value.request('PUT', '/api/v1/works/2611.01070v1/arxiv-id', body={})


def test_unexpected_article_or_version_is_rejected():
    uploader.verify_draft_target({'workIdentifier': uploader.IDENTIFIER, 'versionId': None})
    uploader.verify_draft_target({'workIdentifier': uploader.IDENTIFIER, 'versionId': uploader.EXPECTED_VERSION})
    uploader.verify_draft_target({'workIdentifier': uploader.IDENTIFIER,
                                  'versionId': uploader.PREVIOUS_VERSION,
                                  'status': 'uploading', 'receipt': None})
    for target in ({'workIdentifier': '2611.01071', 'versionId': None},
                   {'workIdentifier': uploader.IDENTIFIER, 'versionId': uploader.PREVIOUS_VERSION},
                   {'workIdentifier': uploader.IDENTIFIER, 'versionId': '2611.01070v3'}):
        with pytest.raises(ValueError):
            uploader.verify_draft_target(target)


def test_receipt_must_identify_existing_article_version_two():
    valid = {'identifier': uploader.IDENTIFIER, 'versionId': uploader.EXPECTED_VERSION,
             'receipt': {'identifier': uploader.IDENTIFIER, 'versionId': uploader.EXPECTED_VERSION,
                         'title': uploader.TITLE}}
    assert uploader.validate_receipt(valid) == valid['receipt']
    assert uploader.validate_receipt(dict(valid, identifier='hexagon:'+uploader.IDENTIFIER)) == valid['receipt']
    for where, key, bad in ((None, 'identifier', '2611.01071'),
                            ('receipt', 'identifier', '2611.01071'),
                            (None, 'versionId', uploader.PREVIOUS_VERSION),
                            ('receipt', 'versionId', uploader.PREVIOUS_VERSION)):
        changed = json.loads(json.dumps(valid))
        (changed if where is None else changed[where])[key] = bad
        with pytest.raises(ValueError):
            uploader.validate_receipt(changed)


def test_remote_source_tree_must_match_every_member():
    manifest = {'source_archive_members': [{'path': 'main.tex', 'sha256': 'abc', 'bytes': 123},
                                          {'path': 'anc/example.py', 'sha256': 'def', 'bytes': 12}]}
    files = [{'path': 'main.tex', 'sha256': 'abc', 'sizeBytes': 123},
             {'path': 'anc/example.py', 'sha256': 'def', 'sizeBytes': 12}]
    uploader.verify_remote_files(files, manifest)
    for changed in (files[:1], files+files[:1], [files[0], dict(files[1], sha256='corrupt')]):
        with pytest.raises(ValueError):
            uploader.verify_remote_files(changed, manifest)


def test_inherited_archive_is_removed_only_from_the_new_revision(monkeypatch, tmp_path):
    state = {'identifier': uploader.IDENTIFIER, 'expected_version': uploader.EXPECTED_VERSION,
             'draft_id': 'new-revision', 'created_by_revision_of_original_deposit': True,
             'committed': False}
    files = [{'id': 'expanded-main', 'kind': 'source', 'relativePath': 'main.tex',
              'sha256': 'main-hash', 'sizeBytes': 10},
             {'id': 'inherited-archive', 'kind': 'source-archive', 'relativePath': 'source.zip',
              'sha256': 'archive-hash', 'sizeBytes': 20}]
    monkeypatch.setattr(uploader, 'PRIVATE', tmp_path)
    monkeypatch.setattr(uploader, 'local_release', lambda: ({}, {}))
    monkeypatch.setattr(uploader, 'state_for', lambda manifest: state)
    monkeypatch.setattr(uploader, 'previous_target', lambda: 'published-original')
    monkeypatch.setattr(uploader, 'save_state', lambda value: None)

    class CloneSession:
        def __init__(self):
            self.remaining = list(files)
            self.deleted = []

        def request(self, method, path):
            if method == 'DELETE':
                assert path == '/api/v1/submissions/new-revision/uploads/inherited-archive'
                self.deleted.append(path)
                self.remaining = []  # The service removes the archive and its members together.
                return {}
            if path == '/api/v1/submissions/published-original':
                return {'status': 'published', 'workIdentifier': uploader.IDENTIFIER,
                        'versionId': uploader.PREVIOUS_VERSION, 'uploads': files}
            assert path == '/api/v1/submissions/new-revision'
            return {'status': 'uploading', 'workIdentifier': uploader.IDENTIFIER,
                    'versionId': uploader.PREVIOUS_VERSION, 'receipt': None,
                    'uploads': self.remaining}

    value = CloneSession()
    uploader.clear_inherited_tree(value)
    assert len(value.deleted) == 1 and state['inherited_v1_file_copy_count'] == 2
    # A state accidentally naming the published parent is rejected before any request.
    state['draft_id'] = 'published-original'
    with pytest.raises(ValueError, match='Only inherited files'):
        uploader.clear_inherited_tree(value)
    assert len(value.deleted) == 1

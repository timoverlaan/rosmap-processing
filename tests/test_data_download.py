"""Tests for pinned, MD5-checked Synapse downloads (no network: Synapse is faked)."""

import hashlib
from pathlib import Path

import pytest

from rosmap_processing.data import download
from rosmap_processing.utils.constants import SynapseFile

CONTENT = b"pinned file content\n"


def _spec(name="a.h5Seurat", content=CONTENT, version=2):
    return SynapseFile(name, version, hashlib.md5(content).hexdigest(), len(content))


@pytest.fixture
def pinned(monkeypatch):
    files = {"syn1": _spec()}
    monkeypatch.setattr(download, "SYNAPSE_FILES", files)
    monkeypatch.setattr(download, "SYNAPSE_OPTIONAL", set())
    return files


class FakeSynapse:
    def __init__(self, content=CONTENT):
        self.content = content
        self.calls = []

    def get(self, syn_id, version, downloadLocation, ifcollision):
        self.calls.append((syn_id, version))
        path = Path(downloadLocation) / "a.h5Seurat"
        path.write_bytes(self.content)
        return type("Entity", (), {"path": str(path), "name": path.name})()


def test_check_file(tmp_path):
    spec = _spec()
    path = tmp_path / spec.name
    assert download.check_file(path, spec) == (False, "missing")
    path.write_bytes(CONTENT)
    assert download.check_file(path, spec) == (True, "")
    path.write_bytes(CONTENT.upper())  # same size, different content
    ok, reason = download.check_file(path, spec)
    assert not ok and reason.startswith("MD5")


def test_downloads_pinned_version_and_verifies(tmp_path, pinned, monkeypatch):
    fake = FakeSynapse()
    monkeypatch.setattr(download, "_login", lambda token: fake)
    download.download_from_synapse({"a": "syn1"}, tmp_path)
    assert fake.calls == [("syn1", 2)]
    assert (tmp_path / "a.h5Seurat").read_bytes() == CONTENT


def test_verified_file_is_kept_without_login(tmp_path, pinned, monkeypatch):
    (tmp_path / "a.h5Seurat").write_bytes(CONTENT)
    monkeypatch.setattr(download, "_login", lambda token: pytest.fail("should not log in"))
    download.download_from_synapse({"a": "syn1"}, tmp_path)


def test_modified_file_is_an_error_not_overwritten(tmp_path, pinned, monkeypatch):
    path = tmp_path / "a.h5Seurat"
    path.write_bytes(b"rewritten locally")
    monkeypatch.setattr(download, "_login", lambda token: pytest.fail("should not download"))
    with pytest.raises(ValueError, match="not overwritten"):
        download.download_from_synapse({"a": "syn1"}, tmp_path)
    assert path.read_bytes() == b"rewritten locally"


def test_corrupt_download_is_an_error(tmp_path, pinned, monkeypatch):
    monkeypatch.setattr(download, "_login", lambda token: FakeSynapse(content=b"truncated"))
    with pytest.raises(ValueError, match="does not match"):
        download.download_from_synapse({"a": "syn1"}, tmp_path)


def test_verify_files_reports_problems(tmp_path, pinned):
    assert download.verify_files({"a": "syn1"}, tmp_path) == ["a.h5Seurat: missing"]
    (tmp_path / "a.h5Seurat").write_bytes(CONTENT)
    assert download.verify_files({"a": "syn1"}, tmp_path) == []


def test_every_configured_file_is_pinned():
    from rosmap_processing.utils.constants import (
        SYNAPSE_FILES, SYNAPSE_IDS_ROSMAP, SYNAPSE_IDS_ROSMAP_MIT,
    )
    for syn_id in [*SYNAPSE_IDS_ROSMAP.values(), *SYNAPSE_IDS_ROSMAP_MIT.values()]:
        assert syn_id in SYNAPSE_FILES

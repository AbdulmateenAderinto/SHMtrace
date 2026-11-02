from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def examples() -> Path:
    return ROOT / "examples"


@pytest.fixture(scope="session")
def wt_fastq(examples) -> Path:
    return examples / "WT3_JH4_sub500.fastq"


@pytest.fixture(scope="session")
def ko_fastq(examples) -> Path:
    return examples / "KO3_JH4_sub500.fastq"


@pytest.fixture(scope="session")
def reference(examples) -> Path:
    return examples / "JH4_reference.fasta"


@pytest.fixture(scope="session")
def contaminant(examples) -> Path:
    return examples / "JH4_contaminant.fasta"


@pytest.fixture(autouse=True)
def _saved_projects_in_tmp(tmp_path, monkeypatch):
    """Runs the web app finishes are saved as projects (session.py); keep them out of the real home."""
    monkeypatch.setenv("SHMTRACE_HOME", str(tmp_path / "shmtrace_home"))

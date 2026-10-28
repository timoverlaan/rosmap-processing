"""Record which code produced an output.

Inside the container there is no .git (only src/ is bind-mounted), so the SLURM
scripts read the commit on the host and pass it in as GIT_COMMIT / GIT_BRANCH /
GIT_DIRTY (see slurm/common.sh). Outside a job, git is asked directly.
"""

import os
import subprocess
from datetime import datetime
from typing import Any, Dict, Optional


def _git(*args: str) -> Optional[str]:
    try:
        return subprocess.run(
            ["git", *args], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_git_info() -> Dict[str, Any]:
    """Commit, branch and number of uncommitted changes, or 'unknown'."""
    commit = os.environ.get("GIT_COMMIT")
    if commit:
        return {
            "commit": commit,
            "branch": os.environ.get("GIT_BRANCH", "unknown"),
            "uncommitted_changes": int(os.environ.get("GIT_DIRTY", "0") or 0),
            "source": "environment",
        }
    commit = _git("rev-parse", "HEAD")
    if commit is None:
        return {"commit": "unknown", "branch": "unknown", "uncommitted_changes": -1, "source": "unavailable"}
    status = _git("status", "--porcelain") or ""
    return {
        "commit": commit,
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD") or "unknown",
        "uncommitted_changes": len([line for line in status.splitlines() if line.strip()]),
        "source": "git",
    }


def provenance_record(step: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """A dict for `adata.uns` describing the step, its parameters and the code version.

    Values are converted to strings/numbers so the record can be written to h5ad.
    """
    from rosmap_processing import __version__

    def _plain(value):
        if value is None:
            return "None"
        if isinstance(value, (bool, int, float, str)):
            return value
        return str(value)

    git = get_git_info()
    return {
        "step": step,
        "package_version": __version__,
        "git_commit": git["commit"],
        "git_branch": git["branch"],
        "git_uncommitted_changes": git["uncommitted_changes"],
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "none"),
        "created": datetime.now().isoformat(timespec="seconds"),
        "params": {k: _plain(v) for k, v in params.items()},
    }

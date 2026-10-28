"""Download ROSMAP data from Synapse, pinned to versions and checked by MD5.

Every file is listed in `utils.constants.SYNAPSE_FILES` with its version, MD5
and size. A file that is already present is checked and kept if it matches; if
it does not match, the download stops with an error rather than silently
keeping or replacing it. Each downloaded file is checked the same way.
"""

import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rosmap_processing.utils.constants import (
    SYNAPSE_FILES,
    SYNAPSE_IDS_ROSMAP,
    SYNAPSE_IDS_ROSMAP_MIT,
    SYNAPSE_OPTIONAL,
    SynapseFile,
)
from rosmap_processing.utils.logging import get_logger

logger = get_logger(__name__)

DATASETS = {
    "rosmap": ("ROSMAP", SYNAPSE_IDS_ROSMAP),
    "rosmap-mit": ("ROSMAP_MIT", SYNAPSE_IDS_ROSMAP_MIT),
}


def md5sum(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def check_file(path: Path, spec: SynapseFile) -> Tuple[bool, str]:
    """(True, "") if the file has the pinned size and MD5, else (False, reason)."""
    if not path.exists():
        return False, "missing"
    size = path.stat().st_size
    if size != spec.size:
        return False, f"size {size} != {spec.size}"
    md5 = md5sum(path)
    if md5 != spec.md5:
        return False, f"MD5 {md5} != {spec.md5}"
    return True, ""


def _wanted(synapse_ids: Dict[str, str], include_optional: bool) -> Dict[str, str]:
    return {
        key: syn_id for key, syn_id in synapse_ids.items()
        if include_optional or syn_id not in SYNAPSE_OPTIONAL
    }


def verify_files(
    synapse_ids: Dict[str, str], directory: Union[str, Path], include_optional: bool = False
) -> List[str]:
    """Check the pinned files in a directory; returns a list of problems (empty if all OK)."""
    directory = Path(directory)
    problems = []
    for key, syn_id in _wanted(synapse_ids, include_optional).items():
        spec = SYNAPSE_FILES[syn_id]
        ok, reason = check_file(directory / spec.name, spec)
        if ok:
            logger.info(f"OK       {spec.name} ({syn_id} v{spec.version})")
        else:
            logger.error(f"MISMATCH {spec.name} ({syn_id} v{spec.version}): {reason}")
            problems.append(f"{spec.name}: {reason}")
    return problems


def download_from_synapse(
    synapse_ids: Dict[str, str],
    output_dir: Path,
    token: Optional[str] = None,
    include_optional: bool = False,
) -> List[Path]:
    """
    Download pinned file versions from Synapse into output_dir.

    Files already present with the pinned MD5 are kept. A present file with a
    different size or MD5 raises ValueError (it is never overwritten).

    Returns
    -------
    List[Path]
        Paths of all pinned files in output_dir
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    wanted = _wanted(synapse_ids, include_optional)
    logger.info(f"{len(wanted)} files for {output_dir}")

    syn = None
    paths = []
    for key, syn_id in wanted.items():
        spec = SYNAPSE_FILES.get(syn_id)
        if spec is None:
            raise ValueError(f"{key} ({syn_id}) has no pinned version in SYNAPSE_FILES")
        dest = output_dir / spec.name
        paths.append(dest)

        if dest.exists():
            ok, reason = check_file(dest, spec)
            if ok:
                logger.info(f"Present and verified: {spec.name} ({syn_id} v{spec.version})")
                continue
            raise ValueError(
                f"{dest} exists but is not {syn_id} v{spec.version} ({reason}). "
                f"It is not overwritten: move it away or download into another directory."
            )

        if syn is None:
            syn = _login(token)
        logger.info(f"Downloading {spec.name} ({syn_id} v{spec.version}, {spec.size / 1e9:.2f} GB)")
        entity = syn.get(syn_id, version=spec.version, downloadLocation=str(output_dir), ifcollision="keep.local")
        got = Path(entity.path) if entity.path else output_dir / entity.name
        if got != dest:
            got.rename(dest)
        ok, reason = check_file(dest, spec)
        if not ok:
            raise ValueError(f"Downloaded {dest} does not match {syn_id} v{spec.version}: {reason}")
        logger.info(f"  -> verified {dest.name}")

    logger.info(f"All {len(paths)} files present and verified in {output_dir}")
    return paths


def _login(token: Optional[str]):
    import synapseclient

    try:
        if token:
            logger.info("Logging in to Synapse with provided token")
            return synapseclient.login(silent=True, authToken=token)
        logger.info("Logging in to Synapse with cached credentials")
        return synapseclient.login(silent=True)
    except Exception as e:
        logger.error(f"Failed to login to Synapse: {e}")
        raise ValueError(
            "Synapse authentication failed. Please provide a valid token "
            "via SYNAPSE_AUTH_TOKEN environment variable or token.txt file."
        )


def _token(synapse_config=None) -> Optional[str]:
    from rosmap_processing.utils.config import SynapseConfig

    return (synapse_config or SynapseConfig()).get_token()


def download_dataset(
    dataset: str,
    output_base: Union[str, Path] = "data/raw",
    synapse_config=None,
    include_optional: bool = False,
) -> List[Path]:
    """Download one dataset ('rosmap' or 'rosmap-mit') into output_base/<ROSMAP|ROSMAP_MIT>."""
    subdir, synapse_ids = DATASETS[dataset]
    logger.info("=" * 60)
    logger.info(f"Downloading {dataset} into {Path(output_base) / subdir}")
    logger.info("=" * 60)

    return download_from_synapse(
        synapse_ids=synapse_ids,
        output_dir=Path(output_base) / subdir,
        token=_token(synapse_config),
        include_optional=include_optional,
    )


def download_rosmap_data(output_dir: Union[str, Path] = "data/raw/ROSMAP", synapse_config=None) -> List[Path]:
    """Download the ROSMAP (Green 2024) files into output_dir."""
    return download_from_synapse(SYNAPSE_IDS_ROSMAP, Path(output_dir), _token(synapse_config))


def download_rosmap_mit_data(output_dir: Union[str, Path] = "data/raw/ROSMAP_MIT", synapse_config=None) -> List[Path]:
    """Download the ROSMAP-MIT (Mathys 2023) files into output_dir."""
    return download_from_synapse(SYNAPSE_IDS_ROSMAP_MIT, Path(output_dir), _token(synapse_config))

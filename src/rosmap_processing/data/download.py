"""Download ROSMAP data from Synapse."""

import synapseclient
from pathlib import Path
from typing import Dict, List, Optional
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from rosmap_processing.utils.constants import (
    SYNAPSE_IDS_ROSMAP,
    SYNAPSE_IDS_ROSMAP_MIT,
)
from rosmap_processing.utils.config import SynapseConfig
from rosmap_processing.utils.logging import get_logger

logger = get_logger(__name__)


def download_from_synapse(
    synapse_ids: Dict[str, str],
    output_dir: Path,
    token: Optional[str] = None,
) -> List[Path]:
    """
    Download files from Synapse.
    
    Parameters
    ----------
    synapse_ids : Dict[str, str]
        Dictionary mapping file names to Synapse IDs
    output_dir : Path
        Directory to save downloaded files
    token : str, optional
        Synapse authentication token. If None, attempts to use cached credentials
        
    Returns
    -------
    List[Path]
        List of paths to downloaded files
        
    Raises
    ------
    ValueError
        If authentication fails
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Downloading {len(synapse_ids)} files from Synapse")
    logger.info(f"Output directory: {output_dir}")
    
    # Login to Synapse
    try:
        if token:
            logger.info("Logging in to Synapse with provided token")
            syn = synapseclient.login(silent=True, authToken=token)
        else:
            logger.info("Logging in to Synapse with cached credentials")
            syn = synapseclient.login(silent=True)
    except Exception as e:
        logger.error(f"Failed to login to Synapse: {e}")
        raise ValueError(
            "Synapse authentication failed. Please provide a valid token "
            "via SYNAPSE_AUTH_TOKEN environment variable or token.txt file."
        )
    
    downloaded_files = []
    
    for name, syn_id in synapse_ids.items():
        logger.info(f"Downloading {name} (syn_id: {syn_id})")
        try:
            entity = syn.get(
                syn_id,
                downloadLocation=str(output_dir),
                ifcollision="keep.local"
            )
            downloaded_path = Path(entity.path)
            downloaded_files.append(downloaded_path)
            logger.info(f"  → Saved to: {downloaded_path.name}")
        except Exception as e:
            logger.error(f"  → Failed to download {name}: {e}")
            raise
    
    logger.info(f"Successfully downloaded {len(downloaded_files)} files")
    return downloaded_files


def download_rosmap_data(
    output_dir: Union[str, Path] = "data/raw/ROSMAP",
    synapse_config: Optional[SynapseConfig] = None,
) -> List[Path]:
    """
    Download ROSMAP dataset from Synapse.
    
    Parameters
    ----------
    output_dir : str or Path, default "data/raw/ROSMAP"
        Directory to save downloaded files
    synapse_config : SynapseConfig, optional
        Synapse configuration. If None, uses default config
        
    Returns
    -------
    List[Path]
        List of paths to downloaded files
    """
    logger.info("="*60)
    logger.info("Downloading ROSMAP dataset")
    logger.info("="*60)
    
    if synapse_config is None:
        synapse_config = SynapseConfig()
    
    token = synapse_config.get_token()
    if token is None:
        logger.warning(
            "No Synapse token found. "
            "Set SYNAPSE_AUTH_TOKEN environment variable or create token.txt"
        )
    
    return download_from_synapse(
        synapse_ids=SYNAPSE_IDS_ROSMAP,
        output_dir=Path(output_dir),
        token=token,
    )


def download_rosmap_mit_data(
    output_dir: Union[str, Path] = "data/raw/ROSMAP_MIT",
    synapse_config: Optional[SynapseConfig] = None,
) -> List[Path]:
    """
    Download ROSMAP-MIT dataset from Synapse.
    
    Parameters
    ----------
    output_dir : str or Path, default "data/raw/ROSMAP_MIT"
        Directory to save downloaded files
    synapse_config : SynapseConfig, optional
        Synapse configuration. If None, uses default config
        
    Returns
    -------
    List[Path]
        List of paths to downloaded files
    """
    logger.info("="*60)
    logger.info("Downloading ROSMAP-MIT dataset")
    logger.info("="*60)
    
    if synapse_config is None:
        synapse_config = SynapseConfig()
    
    token = synapse_config.get_token()
    if token is None:
        logger.warning(
            "No Synapse token found. "
            "Set SYNAPSE_AUTH_TOKEN environment variable or create token.txt"
        )
    
    return download_from_synapse(
        synapse_ids=SYNAPSE_IDS_ROSMAP_MIT,
        output_dir=Path(output_dir),
        token=token,
    )


def download_all_data(
    rosmap_dir: Union[str, Path] = "data/raw/ROSMAP",
    rosmap_mit_dir: Union[str, Path] = "data/raw/ROSMAP_MIT",
    synapse_config: Optional[SynapseConfig] = None,
) -> Dict[str, List[Path]]:
    """
    Download both ROSMAP and ROSMAP-MIT datasets.
    
    Parameters
    ----------
    rosmap_dir : str or Path, default "data/raw/ROSMAP"
        Directory for ROSMAP data
    rosmap_mit_dir : str or Path, default "data/raw/ROSMAP_MIT"
        Directory for ROSMAP-MIT data
    synapse_config : SynapseConfig, optional
        Synapse configuration
        
    Returns
    -------
    Dict[str, List[Path]]
        Dictionary with 'rosmap' and 'rosmap_mit' keys containing file paths
    """
    logger.info("="*60)
    logger.info("Downloading all ROSMAP datasets")
    logger.info("="*60)
    
    rosmap_files = download_rosmap_data(rosmap_dir, synapse_config)
    rosmap_mit_files = download_rosmap_mit_data(rosmap_mit_dir, synapse_config)
    
    return {
        "rosmap": rosmap_files,
        "rosmap_mit": rosmap_mit_files,
    }


if __name__ == "__main__":
    import argparse
    from typing import Union
    
    parser = argparse.ArgumentParser(
        description="Download ROSMAP data from Synapse"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["rosmap", "rosmap-mit", "all"],
        default="all",
        help="Which dataset to download"
    )
    parser.add_argument(
        "--rosmap-dir",
        type=str,
        default="data/raw/ROSMAP",
        help="Output directory for ROSMAP data"
    )
    parser.add_argument(
        "--rosmap-mit-dir",
        type=str,
        default="data/raw/ROSMAP_MIT",
        help="Output directory for ROSMAP-MIT data"
    )
    parser.add_argument(
        "--token",
        type=str,
        help="Synapse authentication token (overrides config)"
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    from rosmap_processing.utils.logging import setup_logging
    setup_logging(level=args.log_level)
    
    # Setup Synapse config
    synapse_config = SynapseConfig()
    if args.token:
        logger.info("Using token from command line argument")
        # Create a temporary config with the provided token
        synapse_config.use_env_token = False
        synapse_config.token_file = None
    
    try:
        if args.dataset == "rosmap":
            download_rosmap_data(args.rosmap_dir, synapse_config if not args.token else None)
        elif args.dataset == "rosmap-mit":
            download_rosmap_mit_data(args.rosmap_mit_dir, synapse_config if not args.token else None)
        else:  # all
            download_all_data(
                args.rosmap_dir,
                args.rosmap_mit_dir,
                synapse_config if not args.token else None
            )
        
        logger.info("Download complete!")
        
    except Exception as e:
        logger.error(f"Download failed: {e}", exc_info=True)
        sys.exit(1)

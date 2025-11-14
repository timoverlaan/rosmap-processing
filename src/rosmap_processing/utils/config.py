"""Configuration management for ROSMAP processing pipeline."""

import os
import yaml
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional
from dataclasses import dataclass, field
from datetime import datetime

from .logging import get_logger

logger = get_logger(__name__)


@dataclass
class ProcessingConfig:
    """Configuration for data processing parameters."""
    
    min_genes: int = 200
    min_cells: int = 5
    n_hvgs: int = 2000
    k_neighbors: int = 30
    n_pca_components: int = 50
    target_sum: float = 1e6
    log_transform: bool = True
    individual_pca: bool = False


@dataclass
class PathsConfig:
    """Configuration for data paths."""
    
    raw_data: Path = Path("data/raw")
    processed: Path = Path("data/processed")
    interim: Path = Path("data/interim")
    metadata: Path = Path("data/metadata")
    output_base: Path = Path("output")  # Base output directory
    logs_base: Path = Path("logs")      # Base logs directory
    
    # These will be set dynamically per run
    output: Optional[Path] = None
    logs: Optional[Path] = None
    
    def __post_init__(self):
        """Convert string paths to Path objects."""
        self.raw_data = Path(self.raw_data)
        self.processed = Path(self.processed)
        self.interim = Path(self.interim)
        self.metadata = Path(self.metadata)
        self.output_base = Path(self.output_base)
        self.logs_base = Path(self.logs_base)
        
        # Set output and logs to base if not specified
        if self.output is None:
            self.output = self.output_base
        else:
            self.output = Path(self.output)
            
        if self.logs is None:
            self.logs = self.logs_base
        else:
            self.logs = Path(self.logs)
    
    def setup_run_directories(
        self, 
        run_name: Optional[str] = None,
        job_id: Optional[str] = None,
        timestamp: Optional[str] = None
    ) -> tuple[Path, Path]:
        """
        Create output and log directories for a specific run.
        
        Parameters
        ----------
        run_name : str, optional
            Name of the run (e.g., "rosmap_processing", "seaad_pipeline")
        job_id : str, optional
            SLURM job ID (if running on cluster)
        timestamp : str, optional
            Timestamp string. If None, uses current time.
            
        Returns
        -------
        output_dir : Path
            Path to output directory for this run
        log_dir : Path
            Path to log directory for this run
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Build directory name
        dir_parts = []
        if run_name:
            dir_parts.append(run_name)
        if job_id:
            dir_parts.append(f"job_{job_id}")
        dir_parts.append(timestamp)
        
        dir_name = "_".join(dir_parts)
        
        # Create directories
        output_dir = self.output_base / dir_name
        log_dir = output_dir / "logs"  # Put logs inside output directory
        
        output_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)
        
        # Update current paths
        self.output = output_dir
        self.logs = log_dir
        
        logger.info(f"Created run directories:")
        logger.info(f"  Output: {output_dir}")
        logger.info(f"  Logs: {log_dir}")
        
        return output_dir, log_dir


@dataclass
class SynapseConfig:
    """Configuration for Synapse access."""
    
    use_env_token: bool = True
    token_file: Optional[Path] = Path("token.txt")
    
    def get_token(self) -> Optional[str]:
        """
        Get Synapse authentication token.
        
        Returns
        -------
        str or None
            Authentication token, or None if not found
        """
        if self.use_env_token:
            token = os.environ.get("SYNAPSE_AUTH_TOKEN")
            if token:
                logger.debug("Using Synapse token from environment variable")
                return token
        
        if self.token_file and Path(self.token_file).exists():
            logger.debug(f"Reading Synapse token from {self.token_file}")
            with open(self.token_file, "r") as f:
                return f.read().strip()
        
        logger.warning("No Synapse token found")
        return None


@dataclass
class LoggingConfig:
    """Configuration for logging."""
    
    level: str = "INFO"
    format: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    log_to_file: bool = True


@dataclass
class Config:
    """Main configuration class."""
    
    dataset_name: str = "ROSMAP"
    dataset_type: str = "single-cell-rnaseq"
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    paths: PathsConfig = field(default_factory=PathsConfig)
    synapse: SynapseConfig = field(default_factory=SynapseConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    
    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "Config":
        """
        Load configuration from YAML file.
        
        Parameters
        ----------
        yaml_path : Path
            Path to YAML configuration file
            
        Returns
        -------
        Config
            Configuration object
            
        Raises
        ------
        FileNotFoundError
            If configuration file doesn't exist
        """
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {yaml_path}")
        
        logger.info(f"Loading configuration from {yaml_path}")
        
        with open(yaml_path, "r") as f:
            config_dict = yaml.safe_load(f)
        
        return cls.from_dict(config_dict)
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "Config":
        """
        Create configuration from dictionary.
        
        Parameters
        ----------
        config_dict : dict
            Configuration dictionary
            
        Returns
        -------
        Config
            Configuration object
        """
        # Extract nested configurations
        processing = ProcessingConfig(**config_dict.get("processing", {}))
        paths = PathsConfig(**config_dict.get("paths", {}))
        synapse = SynapseConfig(**config_dict.get("synapse", {}))
        logging_cfg = LoggingConfig(**config_dict.get("logging", {}))
        
        return cls(
            dataset_name=config_dict.get("dataset_name", "ROSMAP"),
            dataset_type=config_dict.get("dataset_type", "single-cell-rnaseq"),
            processing=processing,
            paths=paths,
            synapse=synapse,
            logging=logging_cfg,
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary.
        
        Returns
        -------
        dict
            Configuration as dictionary
        """
        return {
            "dataset_name": self.dataset_name,
            "dataset_type": self.dataset_type,
            "processing": self.processing.__dict__,
            "paths": {k: str(v) for k, v in self.paths.__dict__.items()},
            "synapse": self.synapse.__dict__,
            "logging": self.logging.__dict__,
        }
    
    def save(self, yaml_path: Path) -> None:
        """
        Save configuration to YAML file.
        
        Parameters
        ----------
        yaml_path : Path
            Path to save configuration file
        """
        yaml_path = Path(yaml_path)
        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(yaml_path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Configuration saved to {yaml_path}")
    
    def setup_run(
        self,
        run_name: Optional[str] = None,
        job_id: Optional[str] = None,
        timestamp: Optional[str] = None,
        save_git_info: bool = True,
        copy_config: bool = True,
        config_source: Optional[Path] = None
    ) -> tuple[Path, Path]:
        """
        Set up directories and metadata for a processing run.
        
        Creates timestamped output and log directories, saves git commit hash,
        and copies the configuration file to the output directory.
        
        Parameters
        ----------
        run_name : str, optional
            Name of the run (e.g., "rosmap_processing", "seaad_pipeline")
        job_id : str, optional
            SLURM job ID (if running on cluster)
        timestamp : str, optional
            Timestamp string. If None, uses current time.
        save_git_info : bool, default=True
            Save git commit hash and status to output directory
        copy_config : bool, default=True
            Copy configuration file to output directory
        config_source : Path, optional
            Source config file to copy. If None, saves current config.
            
        Returns
        -------
        output_dir : Path
            Path to output directory for this run
        log_dir : Path
            Path to log directory for this run
        """
        # Create directories
        output_dir, log_dir = self.paths.setup_run_directories(
            run_name=run_name,
            job_id=job_id,
            timestamp=timestamp
        )
        
        # Save git information
        if save_git_info:
            self._save_git_info(output_dir)
        
        # Copy or save configuration
        if copy_config:
            config_dest = output_dir / "config.yaml"
            if config_source and Path(config_source).exists():
                shutil.copy2(config_source, config_dest)
                logger.info(f"Copied config from {config_source} to {config_dest}")
            else:
                self.save(config_dest)
                logger.info(f"Saved current config to {config_dest}")
        
        return output_dir, log_dir
    
    def _save_git_info(self, output_dir: Path) -> None:
        """
        Save git commit hash and status to output directory.
        
        Parameters
        ----------
        output_dir : Path
            Directory to save git information
        """
        git_info_file = output_dir / "git_info.txt"
        
        try:
            # Get git commit hash
            commit_hash = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                check=True
            ).stdout.strip()
            
            # Get git branch
            branch = subprocess.run(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                capture_output=True,
                text=True,
                check=True
            ).stdout.strip()
            
            # Get git status (check for uncommitted changes)
            status = subprocess.run(
                ["git", "status", "--porcelain"],
                capture_output=True,
                text=True,
                check=True
            ).stdout.strip()
            
            # Write to file
            with open(git_info_file, "w") as f:
                f.write(f"Git Commit Hash: {commit_hash}\n")
                f.write(f"Branch: {branch}\n")
                f.write(f"Timestamp: {datetime.now().isoformat()}\n")
                f.write("\n")
                
                if status:
                    f.write("WARNING: Uncommitted changes detected:\n")
                    f.write(status)
                    f.write("\n")
                    logger.warning("Git repository has uncommitted changes!")
                else:
                    f.write("Repository is clean (no uncommitted changes)\n")
            
            logger.info(f"Git info saved to {git_info_file}")
            logger.info(f"Commit: {commit_hash[:8]} on branch '{branch}'")
            
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"Could not save git info: {e}")
            with open(git_info_file, "w") as f:
                f.write(f"Git information not available: {e}\n")


def load_config(config_path: Optional[Path] = None) -> Config:
    """
    Load configuration from file or use defaults.
    
    Parameters
    ----------
    config_path : Path, optional
        Path to configuration file. If None, uses default configuration
        
    Returns
    -------
    Config
        Configuration object
    """
    if config_path:
        return Config.from_yaml(config_path)
    else:
        logger.info("Using default configuration")
        return Config()

"""Configuration management for ROSMAP processing pipeline."""

import os
import yaml
from pathlib import Path
from typing import Any, Dict, Optional
from dataclasses import dataclass, field

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
    output: Path = Path("output")
    logs: Path = Path("logs")
    
    def __post_init__(self):
        """Convert string paths to Path objects."""
        self.raw_data = Path(self.raw_data)
        self.processed = Path(self.processed)
        self.interim = Path(self.interim)
        self.metadata = Path(self.metadata)
        self.output = Path(self.output)
        self.logs = Path(self.logs)


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

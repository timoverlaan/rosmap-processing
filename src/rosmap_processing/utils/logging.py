"""Logging utilities for the ROSMAP processing pipeline."""

import logging
import sys
from pathlib import Path
from typing import Optional


def setup_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
    format_string: Optional[str] = None,
) -> logging.Logger:
    """
    Configure logging for the application.
    
    Parameters
    ----------
    level : str, default "INFO"
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    log_file : Path, optional
        Path to log file. If None, only logs to console
    format_string : str, optional
        Custom format string. If None, uses default format
        
    Returns
    -------
    logging.Logger
        Configured logger instance
        
    Examples
    --------
    >>> from rosmap_processing.utils.logging import setup_logging
    >>> logger = setup_logging(level="DEBUG")
    >>> logger.info("Processing started")
    """
    # Default format
    if format_string is None:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    # Create formatter
    formatter = logging.Formatter(format_string, datefmt="%Y-%m-%d %H:%M:%S")
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    
    # Setup handlers list
    handlers = [console_handler]
    
    # File handler (optional)
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Configure root logger
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        handlers=handlers,
        force=True,  # Override any existing configuration
    )
    
    # Return logger for this package
    logger = logging.getLogger("rosmap_processing")
    logger.setLevel(getattr(logging, level.upper()))
    
    return logger


def get_logger(name: str = "rosmap_processing") -> logging.Logger:
    """
    Get a logger instance.
    
    Parameters
    ----------
    name : str, default "rosmap_processing"
        Logger name
        
    Returns
    -------
    logging.Logger
        Logger instance
        
    Examples
    --------
    >>> from rosmap_processing.utils.logging import get_logger
    >>> logger = get_logger(__name__)
    >>> logger.info("Module loaded")
    """
    return logging.getLogger(name)

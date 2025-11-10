"""Unit tests for the logging utility module."""

import pytest
from pathlib import Path
from rosmap_processing.utils.logging import setup_logging, get_logger


def test_setup_logging_returns_logger():
    """Test that setup_logging returns a logger object."""
    logger = setup_logging(level="INFO")
    assert logger is not None
    assert hasattr(logger, 'info')
    assert hasattr(logger, 'error')
    assert hasattr(logger, 'debug')


def test_get_logger_returns_logger():
    """Test that get_logger returns a logger with correct name."""
    logger = get_logger("test_module")
    assert logger is not None
    assert logger.name == "test_module"


def test_setup_logging_with_file(temp_output_dir):
    """Test logging to file."""
    log_file = temp_output_dir / "test.log"
    logger = setup_logging(level="INFO", log_file=str(log_file))
    
    logger.info("Test message")
    
    # Close all handlers to release the file
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)
    
    assert log_file.exists()
    content = log_file.read_text()
    assert "Test message" in content


@pytest.mark.parametrize("level", ["DEBUG", "INFO", "WARNING", "ERROR"])
def test_setup_logging_levels(level):
    """Test different logging levels."""
    logger = setup_logging(level=level)
    assert logger is not None

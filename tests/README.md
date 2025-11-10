# Test Suite for ROSMAP Processing Pipeline

This directory contains the test suite for the `rosmap-processing` package.

## Test Structure

```
tests/
├── conftest.py              # Pytest fixtures and configuration
├── test_utils_*.py          # Unit tests for utils modules
├── test_data_*.py           # Unit tests for data processing modules
├── test_core_*.py           # Unit tests for core processing modules
└── test_integration.py      # Integration tests with real data
```

## Running Tests

### Run all tests (except integration tests)
```bash
pixi run pytest
```

### Run all tests including integration tests
```bash
pixi run pytest -m ""
```

### Run only integration tests
```bash
pixi run pytest -m integration
```

### Run tests with coverage
```bash
pixi run pytest --cov=src --cov-report=html
```

### Run specific test file
```bash
pixi run pytest tests/test_utils_logging.py
```

### Run specific test
```bash
pixi run pytest tests/test_utils_logging.py::test_setup_logging_returns_logger
```

## Test Categories

### Unit Tests
Test individual modules and functions in isolation:
- `test_utils_*.py` - Utility modules (logging, constants, validation)
- `test_data_*.py` - Data processing modules (category_fix, metadata)
- `test_core_*.py` - Core processing modules (combine, column_mapping)

### Integration Tests
Test full workflows with real data files:
- `test_integration.py` - End-to-end processing with actual ROSMAP data
- These tests are marked with `@pytest.mark.integration`
- Skip them with: `pytest -m "not integration"`

## Test Fixtures

Common fixtures defined in `conftest.py`:

- `test_data_dir` - Path to test data directory
- `mit_h5ad_path` - Path to MIT immune_cells.h5ad
- `mit_metadata_path` - Path to MIT metadata CSV
- `rosmap_h5ad_path` - Path to ROSMAP astrocytes.h5ad
- `rosmap_metadata_path` - Path to ROSMAP clinical CSV
- `temp_output_dir` - Temporary directory for test outputs (auto-cleaned)
- `small_adata` - Small AnnData object (100 cells) for fast tests
- `medium_adata` - Medium AnnData object (1000 cells) for integration tests

## Test Requirements

- Tests use pytest framework
- Install test dependencies: `pixi install` (already includes pytest)
- Some tests require real data files in `data/raw/`
- Integration tests are automatically skipped if data files are missing

## Writing New Tests

### Unit Test Example
```python
def test_my_function(small_adata):
    \"\"\"Test that my function works correctly.\"\"\"
    result = my_function(small_adata)
    assert result is not None
    assert len(result) > 0
```

### Integration Test Example
```python
@pytest.mark.integration
def test_full_pipeline(mit_h5ad_path, temp_output_dir):
    \"\"\"Test full processing pipeline.\"\"\"
    output_file = temp_output_dir / "output.h5ad"
    process_pipeline(mit_h5ad_path, output_file)
    assert output_file.exists()
```

## Continuous Integration

Tests can be run automatically in CI/CD:
```yaml
# Example GitHub Actions workflow
- name: Run tests
  run: |
    pixi run pytest -m "not integration"
```

## Coverage

Generate coverage report:
```bash
pixi run pytest --cov=src --cov-report=html
open htmlcov/index.html  # View coverage report
```

## Debugging Tests

Run tests with verbose output:
```bash
pixi run pytest -vv
```

Run tests and stop at first failure:
```bash
pixi run pytest -x
```

Run tests with Python debugger on failures:
```bash
pixi run pytest --pdb
```

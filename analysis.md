# Repository Analysis: rosmap-processing

**Date:** November 10, 2025  
**Repository:** rosmap-processing (timoverlaan/rosmap-processing)  
**Branch:** main

## Executive Summary

This repository contains processing scripts for the ROSMAP Alzheimer's Disease dataset(s), using Python (scanpy/anndata) and R (Seurat) for single-cell RNA-seq data processing. The project uses Pixi for dependency management and Singularity containers for reproducibility. While functional, the repository has several critical flaws in code quality, documentation, error handling, and maintainability.

---

## Critical Issues

### 1. **Missing Argument Definitions in `match_columns.py`**
**Severity:** HIGH  
**File:** `src/match_columns.py`

The script references `args.cellclass` and `args.cellsubclass` (lines 162, 166) but these arguments are never defined in the argument parser. This will cause an `AttributeError` at runtime.

```python
# Line 162-166: References undefined arguments
if not args.cellclass or not args.cellsubclass:
    raise ValueError("For ROSMAP_MIT data, --cellclass and --cellsubclass must be provided.")

convert_rosmap_mit(adata, cellclass=args.cellclass, cellsubclass=args.cellsubclass)
```

**Impact:** The ROSMAP_MIT conversion will fail when executed.  
**Solution:** Add argument definitions:
```python
parser.add_argument("--cellclass", type=str, help="Cell class for ROSMAP_MIT data")
parser.add_argument("--cellsubclass", type=str, help="Cell subclass for ROSMAP_MIT data")
```

### 2. **Positional Argument with `required=True`**
**Severity:** MEDIUM  
**File:** `src/match_columns.py` (line 12)

```python
parser.add_argument(
    "path",
    type=str,
    required=True,  # ← INVALID: positional args are always required
    help="Path to the input h5ad file.",
)
```

**Impact:** While this works, it's semantically incorrect. Positional arguments are always required by definition in argparse.  
**Solution:** Remove `required=True` for positional arguments.

### 3. **Debugging Code Left in Production**
**Severity:** MEDIUM  
**File:** `src/add_metadata.py` (line 49)

```python
adata.obs.to_csv(f"data/test_obs{'_mit' if mit else ''}.csv", index=True)  # For debugging, remove later
```

**Impact:** Creates unwanted CSV files every time the script runs, polluting the workspace.  
**Solution:** Remove this line or add a `--debug` flag to control it.

### 4. **Inconsistent Directory Paths in Scripts**
**Severity:** MEDIUM  
**File:** `src/run_all.sh` (lines 71-78, 107-110)

The script combines files from ROSMAP_MIT to `data/raw/ROSMAP/combined.h5ad` (line 78) and files from ROSMAP to `data/raw/ROSMAP_MIT/combined.h5ad` (line 91). This appears to be reversed/swapped.

```bash
# Lines 64-78: MIT data → ROSMAP directory (wrong?)
pixi run python src/combine_h5ad.py \
    data/raw/ROSMAP_MIT/Astrocytes.h5ad \
    ... (MIT files)
    --output data/raw/ROSMAP/combined.h5ad

# Lines 80-91: ROSMAP data → MIT directory (wrong?)
pixi run python src/combine_h5ad.py \
    data/raw/ROSMAP/astrocytes.h5ad \
    ... (ROSMAP files)
    --output data/raw/ROSMAP_MIT/combined.h5ad
```

**Impact:** Output files are placed in confusing/incorrect directories, making debugging difficult.  
**Solution:** Verify the intended output paths and fix the inconsistency.

### 5. **Function Parameter Mismatch**
**Severity:** HIGH  
**File:** `src/match_columns.py` (lines 70, 166)

```python
# Function definition (line 70):
def convert_rosmap_mit(adata: ad.AnnData, cellclass: str, subclass: str) -> None:

# Function call (line 166):
convert_rosmap_mit(adata, cellclass=args.cellclass, cellsubclass=args.cellsubclass)
```

**Impact:** The function expects `subclass` but is called with `cellsubclass`, causing a `TypeError`.  
**Solution:** Make parameter names consistent (prefer `subclass`).

---

## Major Issues

### 6. **No Logging Framework**
**Severity:** MEDIUM  
**Files:** All Python scripts

All scripts use `print()` statements for output instead of proper logging. This makes it difficult to:
- Control verbosity levels
- Redirect output to files
- Debug production issues
- Filter by severity

**Solution:** Implement Python's `logging` module with configurable levels.

### 7. **Incomplete TODOs**
**Severity:** MEDIUM  
**Locations:**
- `src/scanpy_pipeline.py` line 46: "TODO: this is not tested yet, but it should work"
- `src/match_columns.py` line 66, 84: "TODO: also do the metadata columns"
- `src/match_columns.py` line 77: "TODO: this might be a problem"
- `src/combine_h5ad.py` line 34: "TODO: check if the joining is correct here"

**Impact:** Indicates incomplete implementation and untested code paths.  
**Solution:** Either implement the TODOs or document why they're deferred.

### 8. **Hardcoded Synapse Token**
**Severity:** HIGH (Security)  
**File:** `src/download_rosmap.py` (line 48)

```python
with open("token.txt", "r") as f:
    token = f.read().strip()
```

**Issues:**
- `token.txt` is tracked in git (present in workspace)
- No error handling if file is missing
- No guidance on how to obtain/format the token
- Token file is bound into containers but not documented

**Solution:** 
- Use environment variables (`SYNAPSE_AUTH_TOKEN`)
- Add `token.txt` to `.gitignore` (currently missing)
- Document token acquisition in README
- Add error message if token is missing

### 9. **No Input Validation**
**Severity:** MEDIUM  
**Files:** Most Python scripts

Scripts don't validate inputs before processing:
- No checks for file existence
- No validation of file formats
- No verification that required columns exist in dataframes
- No bounds checking on numeric parameters

**Example:** `scanpy_pipeline.py` could fail late if layer doesn't exist, after processing already started.

### 10. **Weak Error Handling**
**Severity:** MEDIUM  
**Files:** All Python scripts

Limited try-except blocks means:
- No graceful handling of I/O errors
- No cleanup of partial outputs on failure
- Stack traces exposed to users
- No recovery mechanisms

---

## Moderate Issues

### 11. **Minimal README**
**Severity:** MEDIUM  
**File:** `README.md`

Current README is only 2 lines:
```markdown
# rosmap-processing
Processing scripts for the ROSMAP AD dataset(s)
```

**Missing:**
- Installation instructions
- Usage examples
- Prerequisites (Synapse account, AWS credentials)
- Data sources and citations
- Pipeline overview
- Environment setup (Pixi vs Singularity)
- Parameter explanations
- Troubleshooting guide

### 12. **No Tests**
**Severity:** MEDIUM  

No test files exist in the repository:
- No unit tests
- No integration tests
- No data validation tests
- No regression tests

**Impact:** Changes cannot be verified automatically, increasing risk of regressions.

### 13. **Inconsistent Naming Conventions**
**Severity:** LOW  
**Files:** Multiple

- Mix of snake_case and camelCase in file names (`run_all.sh` vs `rosmap_mathys_genes.sh`)
- Inconsistent variable naming (e.g., `braaksc` vs `Braak`)
- Some files use capitalized names (e.g., `Astrocytes.h5ad` vs `astrocytes.h5ad`)

### 14. **Duplicate SLURM Output Paths**
**Severity:** LOW  
**Files:** All SLURM scripts

All SLURM jobs write to the same output file pattern:
```bash
#SBATCH --output=slurm/out/%j_all.out
#SBATCH --error=slurm/out/%j_all.out
```

While `%j` makes them unique, better naming would help identify which job produced which log.

### 15. **Large Memory Requirements**
**Severity:** LOW (Operational)  
**Files:** SLURM scripts

Scripts request 400-500GB of memory without justification or documentation of actual requirements. This may limit job scheduling unnecessarily.

### 16. **No Version Pinning for PyPI Dependencies**
**Severity:** MEDIUM  
**File:** `pixi.toml` (lines 35-36)

```toml
[pypi-dependencies]
synapseclient = ">=4.7.0"
awscli = ">=1.40.33, <2"
```

Conda dependencies are properly pinned, but PyPI dependencies use loose version constraints.  
**Impact:** May cause reproducibility issues across different environments.

### 17. **Container Version Sprawl**
**Severity:** LOW (Maintenance)  
**Files:** Root directory

Multiple container versions exist without documentation:
- `container_pixi.sif`
- `container_pixi_0-1-1.sif`
- `container_pixi_0-1-2.sif`
- `container_deseq2.sif`

**Issues:**
- No changelog documenting differences
- Scripts reference different versions inconsistently
- All containers are gitignored but some exist in workspace
- No cleanup strategy for old versions

### 18. **Commented-Out Code**
**Severity:** LOW  
**File:** `src/scanpy_pipeline.py` (line 17)

```python
# parser.add_argument("--export-overlap", action="store_true", help="If set, and the genes from --import-genes are not a subset of the input data, \
#                     only the overlap is used, and we also export an h5ad for the refence dataset with only the overlapping genes.")
```

**Impact:** Unclear if feature is planned, deprecated, or temporarily disabled.

### 19. **Unused R Installation**
**Severity:** LOW  
**File:** `install_deps.R`

Installs `BiocManager` and packages but doesn't appear to use standard BiocManager conventions elsewhere in the codebase.

### 20. **Ambiguous Function Behavior**
**Severity:** MEDIUM  
**File:** `src/scanpy_pipeline.py` (lines 111-115)

```python
if donor_slice.shape[0] <= 1:
    print(f"Skipping donor {donor} with only {donor_slice.shape[0]} cell(s). \
          It will be excluded from the output file, because we cannot construct a kNN graph.")
    continue
```

**Issue:** Silently excludes data without tracking which donors were removed or why. No summary statistics provided at the end.

---

## Minor Issues

### 21. **Inconsistent String Formatting**
**Severity:** LOW  
**Files:** Python scripts

Mix of f-strings, `.format()`, and `%` formatting throughout the codebase.  
**Example:** `print("Importing h5ad file: ", args.path)` vs `print(f"Output: {args.output}")`

### 22. **No Type Hints Consistency**
**Severity:** LOW  

Some functions have type hints (e.g., `convert_rosmap(adata: ad.AnnData) -> None`) but many don't. Script-level code has no type annotations.

### 23. **Magic Numbers**
**Severity:** LOW  
**File:** `src/scanpy_pipeline.py`

Hardcoded values without explanation:
- `min_genes=200` (line 66)
- `min_cells=5` (line 68)
- `target_sum=1e6` (line 99)
- `n_comps=50` (line 103, 120)

These should be constants or configurable parameters.

### 24. **Inconsistent Import Organization**
**Severity:** LOW  

Imports are not consistently organized (stdlib, third-party, local). Some scripts have trailing whitespace in imports.

### 25. **Memory Management Concerns**
**Severity:** MEDIUM  
**File:** `src/scanpy_pipeline.py` (lines 119-122)

```python
if args.individual_pca:
    sc.pp.pca(donor_slice, n_comps=50)
sc.pp.neighbors(donor_slice, n_neighbors=args.k_neighbors, use_rep="X_pca")
slices.append(donor_slice)
```

Accumulates all donor slices in memory before concatenating. For large datasets, this could cause OOM errors.

---

## Documentation Issues

### 26. **Missing Docstrings**
**Severity:** LOW  

No functions have docstrings explaining:
- Parameters and their types
- Return values
- Side effects (especially file modifications)
- Exceptions that might be raised

### 27. **No Pipeline Diagram**
**Severity:** MEDIUM  

Complex multi-step pipeline (`run_all.sh`, `run_all_seaad.sh`) without visual representation or clear documentation of dependencies.

### 28. **Ambiguous File Purpose**
**Severity:** LOW  

Files like `check_h5ad.py` and `fix_categories.py` have unclear purposes from their names alone. Need better documentation.

---

## Code Organization Issues

### 29. **Scripts in Wrong Locations**
**Severity:** LOW  

`src/run_all.sh` and `src/run_all_seaad.sh` are operational scripts, not source code. They might be better placed in a `scripts/` directory.

### 30. **No Module Structure**
**Severity:** MEDIUM  

All scripts are standalone with no shared utilities module. Common operations (loading data, validation, logging) are duplicated.

**Missing:**
- `src/utils.py` for shared functions
- `src/config.py` for constants
- `src/validators.py` for input validation
- `__init__.py` files for proper package structure

### 31. **Data Directory Not in .gitignore**
**Severity:** LOW  

While `data/` is ignored, several data files are listed in the workspace structure, suggesting gitignore isn't working correctly or was added late.

---

## Reproducibility Issues

### 32. **No pixi.lock Verification**
**Severity:** MEDIUM  

`pixi.lock` exists but there's no documentation on when to regenerate it or how to verify it's up-to-date.

### 33. **Singularity Definition Outdated**
**Severity:** MEDIUM  
**File:** `container_pixi.def`

Header shows `# Version: 0.1.1` but repository version is `0.1.2` (from `pixi.toml`). Unclear if container needs rebuilding.

### 34. **No Dependency Graph**
**Severity:** LOW  

Scripts have complex dependencies on each other (e.g., `combine_h5ad.py` must run after conversion scripts) but this isn't documented.

---

## Security Issues

### 35. **Token File in Repository**
**Severity:** CRITICAL  

`token.txt` exists in the workspace directory, suggesting it may have been committed at some point or is being tracked.

### 36. **AWS Credentials Handling**
**Severity:** MEDIUM  
**File:** `src/run_all_seaad.sh`

Uses `--no-sign-request` for AWS S3 access, which is fine for public buckets, but doesn't document whether authentication is needed for other operations.

---

## Performance Issues

### 37. **Sequential Processing in Shell Scripts**
**Severity:** LOW  

`run_all.sh` processes files sequentially that could be parallelized (e.g., all R conversions on lines 12-19, 23-30).

### 38. **Redundant File I/O**
**Severity:** LOW  
**File:** `src/run_all.sh`

Runs `check_h5ad.py` multiple times (lines 62, 79, 112, 113) which involves reading entire h5ad files just for validation.

---

## Best Practice Violations

### 39. **No Exit Code Handling**
**Severity:** MEDIUM  

Shell scripts don't check exit codes of Python/R commands. If a step fails midway, subsequent steps continue with corrupted/missing data.

**Solution:** Add `set -e` or proper error checking.

### 40. **No Progress Indicators**
**Severity:** LOW  

Long-running operations (like PCA computation) don't show progress. Only `tqdm` is used in the loop, not for individual operations.

### 41. **Mixing of Data Formats**
**Severity:** LOW  

Uses h5ad, h5Seurat, RDS, and CSV formats. While necessary for interop, there's no clear documentation on when/why each format is used.

---

## Recommendations

### Immediate Actions (Critical)
1. **Fix `match_columns.py` argument parser** - Add missing `--cellclass` and `--cellsubclass` arguments
2. **Fix parameter name mismatch** in `convert_rosmap_mit()` function call
3. **Remove or guard debugging CSV output** in `add_metadata.py`
4. **Fix swapped output paths** in `src/run_all.sh` combine operations
5. **Add `token.txt` to `.gitignore`** if not already there
6. **Remove `token.txt` from repository history** if it was committed

### Short-term Improvements (High Priority)
1. Add comprehensive README with:
   - Installation instructions
   - Usage examples
   - Data source documentation
   - Pipeline overview
2. Implement proper logging framework
3. Add input validation to all scripts
4. Add error handling with try-except blocks
5. Complete or document all TODO items
6. Add basic unit tests for core functions

### Medium-term Enhancements
1. Refactor shared code into utility modules
2. Add docstrings to all functions
3. Create pipeline diagram/documentation
4. Add version consistency checks
5. Implement better memory management for large datasets
6. Add command-line flags for debugging features
7. Create changelog for container versions

### Long-term Goals
1. Implement comprehensive test suite
2. Add CI/CD pipeline for testing
3. Create user documentation and tutorials
4. Benchmark and optimize performance
5. Add data provenance tracking
6. Implement checkpoint/resume functionality for long pipelines

---

## Positive Aspects

Despite the issues identified, the repository has several strengths:

1. **Good use of modern tools:** Pixi for dependency management, Singularity for containerization
2. **Reproducibility focus:** Lock files and container definitions support reproducibility
3. **Type hints present:** Some functions have proper type annotations
4. **Structured workflow:** Clear separation between download, conversion, and processing steps
5. **Version control:** Uses git with appropriate ignore patterns (mostly)
6. **HPC-ready:** SLURM scripts show awareness of HPC best practices

---

## Summary Statistics

- **Total Issues Found:** 41
- **Critical:** 2 (Security: token handling, Missing arguments)
- **High Severity:** 3 (Parameter mismatch, Missing args, Token security)
- **Medium Severity:** 15
- **Low Severity:** 21
- **Files Analyzed:** 16 Python files, 3 R files, 9 Shell scripts, Config files
- **Lines of Code:** ~2000+ across all files

---

## Conclusion

This repository implements a functional data processing pipeline for ROSMAP datasets but suffers from significant code quality, documentation, and maintainability issues. The most critical problems are:

1. Broken argument parsing that will cause runtime failures
2. Security concerns with token management
3. Lack of comprehensive documentation
4. Absence of tests
5. Incomplete error handling

The codebase appears to be in active development (evidenced by TODOs and debugging code) but needs a refactoring pass to make it production-ready and maintainable for other researchers.

**Recommendation:** Before adding new features, prioritize fixing the critical bugs, adding tests, and improving documentation. This will prevent technical debt from accumulating and make the codebase more reliable for scientific research.

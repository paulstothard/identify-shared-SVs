# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Type hints throughout the codebase for improved type safety
- Comprehensive docstrings for all functions
- `pyproject.toml` for modern Python packaging
- Development dependencies configuration (pytest, black, mypy, ruff)
- `--force` flag to skip interactive prompts for automated pipelines
- This CHANGELOG file to track project changes

### Changed

- Improved code documentation and inline comments

## [0.1.0-beta.1] - Initial Release

### Features

- Core functionality to identify shared structural variants between two VCF files
- Support for position overlap and genotype overlap thresholds
- Optional filtering based on opposing homozygotes
- Automatic VCF compression and indexing
- Progress reporting during processing
- Tab-delimited output file for shared variant identifiers
- Comprehensive README with usage examples
- Test script for validation
- Sample input and output files for testing

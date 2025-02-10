# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-02-12
### Changed
- Regression testing adjusted, validates in Jenkins
- default parameters adjusted according to [GBS-5046](https://jira.oicr.on.ca/browse/GBS-5046)
- changed code that generates input list

## [1.3.0] - 2024-06-25
### Added
[GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.2.2] - 2024-04-09
### Changed
- Updated to a reference built using hg38/p12

## [1.2.1] - 2024-04-04
### Changed
- Changed names in vidarrbuild.json

## [1.2.0] - 2024-03-26
### Added
- Added new task makeCSV and headerFormat, and removed readGroupFormat task
### Changed
- Workflow requires an array of fastq files with read-groups as input. A single fastq-pair (or a single fastq file) must also be inputted as an array.
- runDragen task outputs a merged bam file
- Replaced mode parameter with isRNA parameter (false by default)

## [1.1.0] - 2024-03-05
### Changed
- Changes the zippedOut output into a zipped directory. This ensures that extraction creates a new directory instead of tarbombing the working directory.

## [1.0.0] - 2024-03-04
### Added
- Completes lane level alignments using Dragen
- Supports whole transcriptome alignment

Changelog
----
All changes to versions are noted in this file. This project adheres to [Semantic Versioning]. 

## [0.1.3]

### Added
- Seed option for PRNG and reproducibility
- Testing the entries of changes
- Added Travis CI for build support

### Fixed
- SNP now use binomial random function
- Option `-c` for changes had an extra comma
- Comments about change format were inconsistent
- Bug with incorrect reporting of deletion dna
- Bugs for accuracy of `ref_idx` and `alt_idx` calculation

### Changed
- Moved `get_deletion_offset()` out from inner function

[Semantic Versioning]: http://semver.org/
[0.1.3]: https://pypi.python.org/pypi/simulatesv/0.1.3

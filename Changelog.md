# Changelog

- V 1.5.4.0:
   - updated sequence-formats dependency allows more lenient parsing of pileup-data, now also allowing for arbitrary reference alleles (not just ACTGN). This won't affect calling (reads that support an allele that is not in the SNP-file input are treated as before), but will be less disruptive when parsing pileup-input, for example without a bed-file in samtools.
   - improved error output for parsing problems with pileup-format data. Now only a small part of the problematic chunk is output, hopefully easing error interpretation in such cases
   - output a useful error message if the number of samples passed in --sampleNames is inconsistent with the pileup-input
   - `--samplePopName` now accepts multiple pop-names, separated by comma. The number of pop-names must then match the number of samples.
- V 1.5.3.2: fixed a bug in vcf2eigenstrat that would fail on VCFs with missing Quality values.
- V 1.5.3.1: updated to latest GHC pedantic compilation
- V 1.5.3: Upgraded to sequence-formats 1.7.0 introducing an option for plink popName encoding, and improved pileup-Parsing to allow for skip-reference characters
- V 1.5.2: Fixed a bug with --samplePopName having to be entered after -p or -e. Fixed a bug in the sequence-formats dependency.
- V 1.5.1: Added automatic building
- V 1.5.0: Added support for Plink output
- V 1.4.0.4:
    * Fixed eigenstrat-output in pileupCaller to add a dot after the outputprefix before the file extensions.
    * Updated haskell-stack wrapper scripts for EIGENSOFT and ADMIXTOOLS.
    * Moved unmaintained scripts into unmaintained folder.
- V 1.4.0.3: Updated to new sequence-formats version, now including reading of genetic position from eigenstrat files.
- V 1.4.0.1: Improved README, fixed output bug in genoStats.hs
- V 1.4.0: Added single strand mode, and new triallelic treatment.
- V 1.3.1: Bumped dependency on sequence-formats to new sequence-formats-1.4.0, which includes strand-information in pileup data, as well as rsIds in freqSum to output the correct rsId, and an option to parse chromosomes X, Y and MT.
- V 1.3.0: Lots of refactoring. Lots of testing. Removed some features in vcf2eigenstrat and in pileupCaller, including the option in pileupCaller to call without a SNP file.
- V 1.2.4: normaliseBimWithVCF is ready.
- V 1.2.3 : Adapted to newest sequence-formats. Had to change all the chromosome-related code to the newType Chrom datatype. Also started implementing normaliseBimWithVCF.


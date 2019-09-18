V 1.2.3 : Adapted to newest sequence-formats. Had to change all the chromosome-related code to the newType Chrom datatype. Also started implementing normaliseBimWithVCF.

V 1.2.4: normaliseBimWithVCF is ready.

V 1.3.0: Lots of refactoring. Lots of testing. Removed some features in vcf2eigenstrat and in pileupCaller, including
         the option in pileupCaller to call without a SNP file.

V 1.3.1: Bumped dependency on sequence-formats to new sequence-formats-1.4.0, which includes strand-information in pileup data, as well as 
         rsIds in freqSum to output the correct rsId, and an option to parse chromosomes X, Y and MT.

V 1.4.0: Added single strand mode, and new triallelic treatment.

V 1.4.0.1: Improved README, fixed output bug in genoStats.hs

V 1.4.0.3: Updated to new sequence-formats version, now including reading of genetic position from eigenstrat files.

V 1.4.0.4:
* Fixed eigenstrat-output in pileupCaller to add a dot after the outputprefix before the file extensions.
* Updated haskell-stack wrapper scripts for EIGENSOFT and ADMIXTOOLS.
* Moved unmaintained scripts into unmaintained folder.
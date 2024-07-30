![test_state_windows](https://github.com/hgb-bin-proteomics/CandidateSearch/workflows/Test%20for%20Windows/badge.svg)
![test_state_ubuntu](https://github.com/hgb-bin-proteomics/CandidateSearch/workflows/Test%20for%20Ubuntu%2022.04/badge.svg)
![test_state_macos](https://github.com/hgb-bin-proteomics/CandidateSearch/workflows/Test%20for%20macOS/badge.svg)

# CandidateSearch

Proof-of-concept implementation of a search engine that uses [CandidateVectorSearch](https://github.com/hgb-bin-proteomics/CandidateVectorSearch)
to identify the best peptide candidates for a given mass spectrum. *CandidateSearch* is also the computational backend of the non-cleavable crosslink
search in [MS Annika](https://github.com/hgb-bin-proteomics/MSAnnika). *CandidateSearch* creates the vector encodings of peptides and spectra that are
needed for the sparse matrix search of *CandidateVectorSearch*.

*CandidateSearch* can identify peptide candidates from a given mass spectrum without any precursor ion/mass information and no previous knowledge about
potential fixed or variable modifications. *CandidateSearch* can also identify peptidoform candidates if a set of fixed and variable modifications
is provided. The aim of *CandidateSearch* is to reduce the search space for a given identification task by filtering out unlikely peptide or peptidoform
candidates. It is **NOT** meant to be a standalone search engine for peptide/peptidoform identification.

A simplified break down of the *CandidateSearch* algorithm is given in the following:
- Read the given MS2 spectra from the mgf file.
- Generate the encoding vectors for each spectrum.
- Transform spectrum encoding vectors into the representation needed for *CandidateVectorSearch*.
- Read the given fasta file.
- Digest the proteins of the fasta file into peptides.
- [Optional] Generate decoy peptides and peptidoforms.
- Calculate theoretical ion m/z values for all peptides.
- Generate the encoding vectors for each peptide.
- Transform the peptide encoding vectors into the representation needed for *CandidateVectorSearch*.
- Run *CandidateVectorSearch*.
- Process the results of *CandidateVectorSearch*.
- Create a csv file that maps every spectrum (scan number) to a list of the best *n* peptide candidates.
- Done!

## Usage

Running *CandidateSearch* requires three files:
- An mgf file containing MS2 spectra.  
  **We highly recommend to deisotope and deconvolute spectra before search!**
- A fasta file containing sample proteins.
- A settings file containing parameters for digestion, ion calculation and search (see below for an explanation of the settings file).

The *CandidateSearch* executable can then be run like this:

```bash
CandidateSearch.exe spectra.mgf database.fasta settings.txt
```

Example files that can be used to test *CandidateSearch* can be found in `/data`.

## Settings

The settings file accepts the following parameters:
- MAX_CLEAVAGES: The maximum number of allowed missed cleavages during digestion. (integer, default = 2)
- MIN_PEP_LENGTH: The minimum length of a peptide to be considered for search. (integer, default = 5)
- MAX_PEP_LENGTH: The maximum length of a peptide to be considered for search. (integer, default = 30)
- MAX_PRECURSOR_CHARGE: The maximum considered precursor ion charge. (integer, default = 4)
- MAX_FRAGMENT_CHARGE: The maximum considered fragment ion charge. (string, default = +1)
- MAX_NEUTRAL_LOSSES: The maximum number of neutral losses considered during ion calculation. (integer, default = 1)
- MAX_NEUTRAL_LOSS_MODS: The maximum number of neutral loss modifications considered during ion calculation. (integer, default = 2)
- FIXED_MODIFICATIONS: Fixed modifications that should be considered during search given as `(char)amino_acid:(double)modification_mass`. An example
would be carbamidomethylation of cysteine, which would be denoted as `C:57.021464;`. Several fixed modifications can be provided. (string, default = None)
- VARIABLE_MODIFICATIONS: Variable modifications that should be considered during search given as `(char)amino_acid:(double)modification_mass`. An example
would be oxidation of methionine, which would be denoted as `M:15.994915;`. Several variable modifications can be provided. If no modifications are
given, *CandidateSearch* will return the best scoring unmodified peptidoforms for a given spectrum. (string, default = None)
- DECOY_SEARCH: Whether decoy search should be performed or not. Accepts `true` or `false`. (bool, default = true)
- TOP_N: The number of best candidates that should be returned by the search. (integer, default = 1000)
- TOLERANCE: Tolerance used for matching theoretical ions to experimental peaks. Given in Dalton. (double, default = 0.02)
- NORMALIZE: Whether or not *CandidateVectorSearch* scores should be normalized before selecting the best *n* candidates. Accepts `true` or `false`.
(bool, default = false)
- USE_GAUSSIAN: Whether or not experimental peaks should be modelled as gaussian distributions with `mu = (m/z)` and `sigma = (tolerance/3)`.
Accepts `true` or `false`. (bool, default = true)  
- MODE: Search approach used by *CandidateVectorSearch*. One of the following (default = CPU_SMi32):
  - CPU_DVi32: Sparse int matrix - dense int vector search on the CPU.
  - CPU_DVf32: Sparse float matrix - dense float vector search on the CPU.
  - CPU_DMi32: Sparse int matrix - dense int matrix search on the CPU.
  - CPU_DMf32: Sparse float matrix - dense float matrix search on the CPU.
  - CPU_SVi32: Sparse int matrix - sparse int vector search on the CPU.
  - CPU_SVf32: Sparse float matrix - sparse float vector search on the CPU.
  - CPU_SMi32: Sparse int matrix - sparse int matrix search on the CPU.
  - CPU_SMf32: Sparse float matrix - sparse float matrix search on the CPU.
  - GPU_DVf32: Sparse float matrix - dense float vector search on the GPU (see [requirements](https://github.com/hgb-bin-proteomics/CandidateVectorSearch/blob/master/README.md)).
  - GPU_DMf32: Sparse float matrix - dense float matrix search on the GPU (see [requirements](https://github.com/hgb-bin-proteomics/CandidateVectorSearch/blob/master/README.md)).
  - GPU_SMf32: Sparse float matrix - sparse float matrix search on the GPU (see [requirements](https://github.com/hgb-bin-proteomics/CandidateVectorSearch/blob/master/README.md)).

For the last five parameters you might additionally want to check the documentation of
[CandidateVectorSearch](https://github.com/hgb-bin-proteomics/CandidateVectorSearch) to get a better understanding of their meaning.

An empty `settings.txt` file is a valid configuration for search (default parameters will be used), however not providing a settings file at all is
not valid.

An example `settings.txt` file is provided [here](https://github.com/hgb-bin-proteomics/CandidateSearch/blob/master/settings.txt).

Additionally its contents are listed below, which should help in understanding the formatting:

```
## DIGESTIONS PARAMETERS
MAX_CLEAVAGES = 2
MIN_PEP_LENGTH = 5
MAX_PEP_LENGTH = 30

## ION CALCULATION PARAMETERS
MAX_PRECURSOR_CHARGE = 4
MAX_FRAGMENT_CHARGE = +1
MAX_NEUTRAL_LOSSES = 1
MAX_NEUTRAL_LOSS_MODS = 2
#FIXED_MODIFICATIONS = None
FIXED_MODIFICATIONS = C:57.021464;
#VARIABLE_MODIFICATIONS = None
VARIABLE_MODIFICATIONS = M:15.994915;
#VARIABLE_MODIFICATIONS = M:15.994915;K:284.173607;

## SEARCH PARAMETERS
DECOY_SEARCH = true

## VECTOR SEARCH PARAMETERS
TOP_N = 1000
TOLERANCE = 0.02
NORMALIZE = false
USE_GAUSSIAN = true
MODE = CPU_SMi32
```

## Documentation

The code of this search engine is fully documented within the `.cs` code files. A good entry point is the main function of *CandidateSearch* which is
implemented [here](https://github.com/hgb-bin-proteomics/CandidateSearch/blob/master/CandidateSearch.cs). Documentation generated by
[Doxygen](https://github.com/doxygen/doxygen) is also available here:
[https://hgb-bin-proteomics.github.io/CandidateSearch/](https://hgb-bin-proteomics.github.io/CandidateSearch/)

## Requirements

- [.NET](https://dotnet.microsoft.com/en-us/) may be required on Windows systems.
- \[Optional\] Using GPU based approaches requires a CUDA capable GPU and CUDA version == 12.2.0
([download here](https://developer.nvidia.com/cuda-12-2-0-download-archive)). Other CUDA versions may or may not produce the desired results
([see this issue](https://github.com/hgb-bin-proteomics/CandidateVectorSearch/issues/32)).

## Downloads

Compiled DLLs and and executables are available in the `exe` folder or in
[Releases](https://github.com/hgb-bin-proteomics/CandidateSearch/releases).

We supply compiled executables and DLLs for:
- Windows 10/11 (x86, 64-bit)
- Ubuntu 22.04 (x86, 64-bit)
- macOS 14.4 (arm, 64-bit)

For other operating systems/architectures please compile the source code yourself!
You will also need to compile [CandidateVectorSearch](https://github.com/hgb-bin-proteomics/CandidateVectorSearch)!

## Limitations

This a proof-of-concept implementation that shows the applicability of our *CandidateVectorSearch* approach and not a fully fledged search engine,
therefore this implementation comes with a few limitations:

- We currently only have implemented tryptic digestion.
  - You can implement your own digestion [here](https://github.com/hgb-bin-proteomics/CandidateSearch/blob/master/amanda/FASTAParser.cs#L17-L22).
- We currently have not implemented support for N- or C-terminal modifications.
- We currently have only implemented support for one possible modification per amino acid.
- We only support spectra in centroid mode (we can't really do anything with spectra in profile mode).
- We only support databases up to a size of 12 500 000 peptides/peptidoforms, beyond that we can't guarantee that the matrix can be allocated anymore.
  - Consider splitting your fasta into smaller chunks and searching them separately if the generated database size exceeds 12 500 000.
- The limitations of [CandidateVectorSearch](https://github.com/hgb-bin-proteomics/CandidateVectorSearch?tab=readme-ov-file#limitations) also apply here.

## Results

Example results of *CandidateSearch* and results analysis are given in `tests`. An extensive report is given in [results.md](results.md).

![Results on a HeLa dataset](tests/v1.0.0/results.svg)

**Figure 1:** Identifying peptide candidates and peptidoform candidates with *CandidateSearch [v1.0.0]* in a
[HeLa dataset](https://www.ebi.ac.uk/pride/archive/projects/PXD007750) using the human swissprot database. The considered ground truth was
an [MS Amanda](https://ms.imp.ac.at/?goto=msamanda) search validated with [Percolator](https://github.com/percolator/percolator). For every
high-confidence PSM we checked if the identified peptide/peptidoform was among the top 50/100/500/1000 hits of *CandidateSearch*. We reach
almost 100% coverage within the first 1000 hits of *CandidateSearch* (for reference: the whole database contained ~4 200 000 peptides or ~10 500 000 peptidoforms).

## Benchmarks

Benchmarks of the different algorithms can be found in [benchmarks.md](benchmarks.md).

![benchmark_hpc_1A](benchmarks/vis_A/1A.svg)

**Figure 2:** Int32-based sparse matrix * dense matrix search using
[Eigen](https://eigen.tuxfamily.org/) generally yields the fastest computation
time on modern CPUs.

## Known Issues

[List of known issues](https://github.com/hgb-bin-proteomics/CandidateSearch/issues)

## Citing

If you are using [parts of] *CandidateSearch* please cite:

```
MS Annika 3.0 (publication wip)
```

## License

- [MIT](https://github.com/hgb-bin-proteomics/CandidateSearch/blob/master/LICENSE)

## Contact

[micha.birklbauer@fh-hagenberg.at](mailto:micha.birklbauer@fh-hagenberg.at)

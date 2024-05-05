# Benchmarks

The following are benchmarks of the different sparse matrix/vector multiplication
methods of [Eigen](https://eigen.tuxfamily.org/) and
[cuSPARSE](https://docs.nvidia.com/cuda/cusparse/) that are implemented in
_CandidateSearch_ using real mass spectrometry data.

We ran every benchmark five times to get a more comprehensive overview of
computation times. The averages are plotted below, with error bars denoting
standard deviation.

For all benchmarks we search the spectra from file
`benchmarks/XLpeplib_Beveridge_QEx-HFX_DSS_R1_deconvoluted.mgf` against the
database `benchmarks/Cas9+uniprotkb_proteome_UP000005640_AND_revi_2024_03_22.fasta`
(Cas9 + human SwissProt sequences) with settings `benchmarks/settings.txt`. All
benchmarks were conducted during light background usage (e.g. open browser, text
editor, etc.).

## Abbreviations

The following terms are used synonymously throughout the document:
- `f32CPU_SV`: Float32-(CPU-)based sparse matrix * sparse vector search (using [Eigen](https://eigen.tuxfamily.org/))
- `i32CPU_SV`: Int32-(CPU-)based sparse matrix * sparse vector search (using [Eigen](https://eigen.tuxfamily.org/))
- `f32CPU_DV`: Float32-(CPU-)based sparse matrix * dense vector search (using [Eigen](https://eigen.tuxfamily.org/))
- `i32CPU_DV`: Int32-(CPU-)based sparse matrix * dense vector search (using [Eigen](https://eigen.tuxfamily.org/))
- `f32CPU_SM`: Float32-(CPU-)based sparse matrix * sparse matrix search (using [Eigen](https://eigen.tuxfamily.org/))
- `i32CPU_SM`: Int32-(CPU-)based sparse matrix * sparse matrix search (using [Eigen](https://eigen.tuxfamily.org/))
- `f32CPU_DM`: Float32-(CPU-)based sparse matrix * dense matrix search (using [Eigen](https://eigen.tuxfamily.org/))
- `i32CPU_DM`: Int32-(CPU-)based sparse matrix * dense matrix search (using [Eigen](https://eigen.tuxfamily.org/))
- `f32GPU_DV`: Float32-(GPU-)based sparse matrix * dense vector search (using [cuSPARSE](https://docs.nvidia.com/cuda/cusparse/))
- `f32GPU_DM`: Float32-(GPU-)based sparse matrix * dense matrix search (using [cuSPARSE](https://docs.nvidia.com/cuda/cusparse/))
- `f32GPU_SM`: Float32-(GPU-)based sparse matrix * sparse matrix search (using [cuSPARSE](https://docs.nvidia.com/cuda/cusparse/))

## Hardware

The system we tested this on was a desktop PC with the following hardware:
- MB: ASUS ROG Strix B650E-I
- CPU: AMD Ryzen 7900X [12 cores @ 4.7 GHz base / 5.6 GHz boost]
- RAM: Kingston 64 GB DDR5 RAM [5600 MT/s, 36 CAS]
- GPU: ASUS Dual [Nvidia] GeForce RTX 4060 Ti OC [16 GB VRAM]*
- SSD/HDD: Corsair MP600 Pro NH 2 TB NVMe SSD [PCIe 4.0]
- OS: Windows 11 Pro 64-bit (10.0, Build 22631)

*_Note:_ `Dual` _is part of the name, this is a single graphics card!_

### Normalization = False && Use_Gaussian = False

### Normalization = False && Use_Gaussian = True

### Normalization = True && Use_Gaussian = False

### Normalization = True && Use_Gaussian = True

## Conclusions

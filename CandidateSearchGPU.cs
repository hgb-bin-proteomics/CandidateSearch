using CandidateSearch.util;
using System.Diagnostics;

namespace CandidateSearch
{
    /// <summary>
    /// Searching for candidates on the GPU.
    /// </summary>
    public static class CandidateSearchGPU
    {
        /// <summary>
        /// Searches the given spectra for peptide candidates of the given database using the given settings on the GPU.
        /// </summary>
        /// <param name="spectraFile">Filename of the mgf file containing the MS2 spectra.</param>
        /// <param name="databaseFile">Filename of the fasta file containing the proteins that should be considered for search.</param>
        /// <param name="settings">Settings containing parameters for digestion, ion calculation and search.</param>
        public static void Search(string spectraFile, string databaseFile, Settings settings)
        {
            var spectra = MGFReader.readMGF(spectraFile);
            Console.WriteLine($"Read {spectra.Count} spectra.");

            var peptides = DatabaseReader.readFASTA(databaseFile, settings, generateDecoys: settings.DECOY_SEARCH);
            Console.WriteLine($"Generated {peptides.Count} peptides/peptidoforms from fasta file.");

            if (peptides.Count > 12500000)
            {
                Console.WriteLine("Database size exceeds 12 500 000, might not be able to allocate a matrix of that size!");
                Console.WriteLine("Please see 'Limitations' for more info and possible work arounds!");
            }

            // sorting database
            Console.WriteLine("Sorting peptides/peptidoforms...");
            var sortTime = new Stopwatch(); sortTime.Start();
            peptides.Sort((x, y) => x.ToString().CompareTo(y.ToString()));
            sortTime.Stop();
            Console.WriteLine($"Sorted peptides/peptidoforms for search in {sortTime.Elapsed.TotalSeconds} seconds.");

            // generating the csrColIdx and csrRowoffsets arrays
            int csrColIdxLength = 0;
            foreach (var peptide in peptides)
            {
                var encoding = peptide.getEnconding();
                csrColIdxLength += encoding.Length;
            }

            var csrColIdx = new int[csrColIdxLength];
            var csrRowoffsets = new int[peptides.Count + 1];

            int currentIdxCsrColIdx = 0;
            int currentIdxCsrRowoffsets = 0;
            foreach (var peptide in peptides)
            {
                csrRowoffsets[currentIdxCsrRowoffsets] = currentIdxCsrColIdx;
                currentIdxCsrRowoffsets++;
                var encoding = peptide.getEnconding();
                foreach (var value in encoding)
                {
                    csrColIdx[currentIdxCsrColIdx] = value;
                    currentIdxCsrColIdx++;
                }
            }

            csrRowoffsets[peptides.Count] = csrColIdxLength;

            // generating the spectraValues and spectraIdx arrays
            int spectraValuesLength = 0;
            foreach (var spectrum in spectra)
            {
                var encoding = spectrum.getEncoding();
                spectraValuesLength += encoding.Length;
            }

            var spectraValues = new int[spectraValuesLength];
            var spectraIdx = new int[spectra.Count];

            int currentIdxSV = 0;
            int currentIdxSI = 0;
            foreach (var spectrum in spectra)
            {
                spectraIdx[currentIdxSI] = currentIdxSV;
                currentIdxSI++;
                var encoding = spectrum.getEncoding();
                foreach (var value in encoding)
                {
                    spectraValues[currentIdxSV] = value;
                    currentIdxSV++;
                }
            }

            VectorSearchInterface.VectorSearchAPI.GPU_METHODS METHOD;
            switch (settings.MODE)
            {
                case "GPU_DVf32":
                    METHOD = VectorSearchInterface.VectorSearchAPI.GPU_METHODS.f32GPU_DV;
                    break;
                case "GPU_DMf32":
                    METHOD = VectorSearchInterface.VectorSearchAPI.GPU_METHODS.f32GPU_DM;
                    break;
                case "GPU_SMf32":
                    METHOD = VectorSearchInterface.VectorSearchAPI.GPU_METHODS.f32GPU_SM;
                    break;
                default:
                    METHOD = VectorSearchInterface.VectorSearchAPI.GPU_METHODS.f32GPU_DV;
                    break;
            }

            var sw = new Stopwatch();
            sw.Start();

            var result = VectorSearchInterface.VectorSearchAPI.searchGPU(ref csrRowoffsets,
                                                                         ref csrColIdx,
                                                                         ref spectraValues,
                                                                         ref spectraIdx,
                                                                         topN: settings.TOP_N,
                                                                         tolerance: settings.TOLERANCE,
                                                                         normalize: settings.NORMALIZE,
                                                                         useGaussianTol: settings.USE_GAUSSIAN,
                                                                         batchSize: 100,
                                                                         method: METHOD,
                                                                         verbose: 1000,
                                                                         memStat: out int memStat);

            sw.Stop();

            Console.WriteLine($"GPU search finished with code {memStat}. Search took {sw.Elapsed.TotalSeconds} seconds.");

            var processedResult = new Result(ref result, ref peptides, ref spectra, TopN: settings.TOP_N);
            var csvStat = processedResult.export(spectraFile + "_results.csv");

            Console.WriteLine($"Result file written to disk with code {csvStat}. Search finished!");

            return;
        }
    }
}

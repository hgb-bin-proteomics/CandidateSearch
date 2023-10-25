using CandidateSearch.util;
using System.Diagnostics;

namespace CandidateSearch
{
    public static class CandidateSearchGPU
    {
        public static void Search(string spectraFile, string databaseFile, Settings settings)
        {
            var spectra = MGFReader.readMGF(spectraFile, settings.DECONVOLUTE_SPECTRA);
            Console.WriteLine($"Read {spectra.Count} spectra.");

            var peptides = DatabaseReader.readFASTA(databaseFile, settings, generateDecoys: settings.DECOY_SEARCH);
            Console.WriteLine($"Generated {peptides.Count} peptides from fasta file.");

            int candidateValuesLength = 0;
            foreach (var peptide in peptides)
            {
                var encoding = peptide.getEnconding();
                candidateValuesLength += encoding.Length;
            }

            var candidateValues = new int[candidateValuesLength];
            var candidateIdx = new int[peptides.Count + 1];

            int currentIdxCV = 0;
            int currentIdxCI = 0;
            foreach (var peptide in peptides)
            {
                candidateIdx[currentIdxCI] = currentIdxCV;
                currentIdxCI++;
                var encoding = peptide.getEnconding();
                foreach (var value in encoding)
                {
                    candidateValues[currentIdxCV] = value;
                    currentIdxCV++;
                }
            }

            candidateIdx[peptides.Count] = candidateValuesLength;

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

            var sw = new Stopwatch();
            sw.Start();

            var result = VectorSearchInterface.VectorSearchAPI.searchGPU(ref candidateValues,
                                                                         ref candidateIdx,
                                                                         ref spectraValues,
                                                                         ref spectraIdx,
                                                                         topN: settings.TOP_N,
                                                                         tolerance: settings.TOLERANCE,
                                                                         normalize: settings.NORMALIZE,
                                                                         useGaussianTol: settings.USE_GAUSSIAN,
                                                                         batched: settings.MODE == "GPU_DM" || settings.MODE == "GPU_SM",
                                                                         batchSize: 100,
                                                                         useSparse: settings.MODE == "GPU_SM",
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

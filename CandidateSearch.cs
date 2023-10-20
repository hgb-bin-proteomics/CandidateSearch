using CandidateSearch.util;
using System.Diagnostics;

namespace CandidateSearch
{
    public class CandidateSearch
    {
        // config digestion
        public const int MAX_CLEAVAGES = 2;
        public const int MIN_PEP_LENGTH = 5;
        public const int MAX_PEP_LENGTH = 30;

        // config ion calculation
        public const int MAX_CHARGE = 3;
        public const int MAX_NEUTRAL_LOSSES = 1;
        public const int MAX_NEUTRAL_LOSS_MODS = 2;
        public const string MAX_ALLOWED_CHARGE = "+3";

        // config vector search
        public const int TOP_N = 100;
        public const float TOLERANCE = 0.02f;
        public const bool NORMALIZE = true;
        public const bool USE_GAUSSIAN = true;

        public static void Main(string[] args)
        {
            if (args.Length == 3) {
                var spectraFile = args[0];
                var databaseFile = args[1];
                var settingsFile = args[2];

                var spectra = MGFReader.readMGF(spectraFile);
                Console.WriteLine($"Read {spectra.Count} spectra.");

                var peptides = DatabaseReader.readFASTA(databaseFile, generateDecoys: true);
                Console.WriteLine($"Generated {peptides.Count} peptides from fasta file.");

                int candidateValuesLength = 0;
                foreach (var peptide in peptides)
                {
                    var encoding = peptide.getEnconding();
                    candidateValuesLength += encoding.Length;
                }

                var candidateValues = new int[candidateValuesLength];
                var candidateIdx = new int[peptides.Count];

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
                foreach(var spectrum in spectra)
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

                int memStat;
                var result = VectorSearchInterface.VectorSearchAPI.searchCPU(ref candidateValues, 
                                                                             ref candidateIdx, 
                                                                             ref spectraValues, 
                                                                             ref spectraIdx,
                                                                             topN: TOP_N, 
                                                                             tolerance: TOLERANCE, 
                                                                             normalize: NORMALIZE, 
                                                                             useGaussianTol: USE_GAUSSIAN,
                                                                             batched: false, 
                                                                             batchSize: 100, 
                                                                             useSparse: false, 
                                                                             cores: 0, 
                                                                             verbose: 1000,
                                                                             memStat: out memStat);

                sw.Stop();

                Console.WriteLine($"CPU search finished with code {memStat}. Search took {sw.Elapsed.TotalSeconds} seconds.");

                var processedResult = new Result(ref result, ref peptides, ref spectra, TopN: 100);
                var csvStat = processedResult.export(spectraFile + "_results.csv");

                Console.WriteLine($"Result file written to disk with code {csvStat}. Search finished!");

                return;
            }

            Console.WriteLine("Incorrect number of arguments! CandidateSearch needs exactly 3 arguments: spectra.mgf database.fasta settings.xml");
            return;
        }
    }
}
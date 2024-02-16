using CandidateSearch.util;

namespace CandidateSearch
{
    /// <summary>
    /// Searching for candidates using VectorSearch.
    ///     Requires VectorSearch.dll     v1.7.2
    ///     Requires VectorSearchCUDA.dll v1.4.8
    /// </summary>
    public class CandidateSearch
    {
        /// <summary>
        /// Current version of CandidateSearch.
        /// </summary>
        public const string version = "1.1.1";

        /// <summary>
        /// Main function that is executed when CandidateSearch is run.
        /// </summary>
        /// <param name="args">Arguments passed via commandline.</param>
        public static void Main(string[] args)
        {
            if (args.Length == 3) {
                var spectraFile = args[0];
                var databaseFile = args[1];
                var settingsFile = args[2];

                Console.WriteLine($"Starting Candidate Search v{version} ...");

                var settings = SettingsReader.readSettings(settingsFile);
                Console.WriteLine($"Read settings file '{settingsFile}' with the following settings:");
                Console.WriteLine(settings.ToString());

                if (settings.MODE.Split("_").First().Trim() == "GPU")
                {
                    CandidateSearchGPU.Search(spectraFile, databaseFile, settings);
                }
                else
                {
                    CandidateSearchCPU.Search(spectraFile, databaseFile, settings);
                }

                return;
            }

            Console.WriteLine("Incorrect number of arguments! CandidateSearch needs exactly 3 arguments: spectra.mgf database.fasta settings.txt");
            return;
        }
    }
}
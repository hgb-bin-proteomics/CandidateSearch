using CandidateSearch.util;

namespace CandidateSearch
{
    public class CandidateSearch
    {

        public static void Main(string[] args)
        {
            if (args.Length == 3) {
                var spectraFile = args[0];
                var databaseFile = args[1];
                var settingsFile = args[2];

                var spectra = MGFReader.readMGF(spectraFile);
                var peptides = DatabaseReader.readFASTA(databaseFile);

                return;
            }

            Console.WriteLine("Incorrect number of arguments! CandidateSearch needs exactly 3 arguments: spectra.mgf database.fasta settings.xml");
            return;
        }
    }
}
using CandidateSearch.util;

namespace CandidateSearch
{
    public class CandidateSearch
    {

        public static void Main(string[] args)
        {
            if (args.Length == 3) { 
                var spectra = MGFReader.readMGF(args[0]);
                var database = args[1];
                var settings = args[2];

                return;
            }

            Console.WriteLine("Incorrect number of arguments! CandidateSearch needs exactly 3 arguments: spectra.mgf database.fasta settings.xml");
            return;
        }
    }
}
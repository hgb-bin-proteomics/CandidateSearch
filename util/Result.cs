namespace CandidateSearch.util
{
    public class Result
    {
        public Dictionary<int, List<Peptide>> result;

        public Result(ref int[] SearchResult, ref List<Peptide> Peptides, ref List<Spectrum> Spectra, int TopN) {
            
            result = new Dictionary<int, List<Peptide>>();
            
            int currentSearchResultIdx = 0;
            foreach (var spectrum in Spectra)
            {
                result.Add(spectrum.scanNumber, new List<Peptide>());
                for (int i = currentSearchResultIdx; i < currentSearchResultIdx + TopN; i++)
                {
                    result[spectrum.scanNumber].Add(Peptides[SearchResult[i]]);
                }

                currentSearchResultIdx += TopN;
            }
        }

        public int export(string filename)
        {
            int status = 1;

            try
            {
                // code
            }
            catch (Exception ex)
            {
                Console.WriteLine("Something went wrong:");
                Console.WriteLine(ex.ToString());
                status = 1;
            }

            return status;
        }
    }
}

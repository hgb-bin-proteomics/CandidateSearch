namespace CandidateSearch.util
{
    /// <summary>
    /// Result class to transform and store the returned integer array of VectorSearch.
    /// </summary>
    public class Result
    {
        /// <summary>
        /// Dictionary that maps scan numbers to lists of peptide candidates.
        /// </summary>
        public Dictionary<int, List<Peptide>> result { get; }

        /// <summary>
        /// Constructor creating a result item that stores the processed VectorSearch result. 
        /// </summary>
        /// <param name="SearchResult">The search result of the VectorSearch.</param>
        /// <param name="Peptides">The complete list of considered peptides/peptidoforms.</param>
        /// <param name="Spectra">The complete list of searched spectra.</param>
        /// <param name="TopN">The number of top candidates reported for each spectrum.</param>
        public Result(ref int[] SearchResult, ref List<Peptide> Peptides, ref List<Spectrum> Spectra, int TopN) {
            
            result = new Dictionary<int, List<Peptide>>();
            
            int currentSearchResultIdx = 0;
            foreach (var spectrum in Spectra)
            {
                if (!result.ContainsKey(spectrum.scanNumber))
                {
                    result.Add(spectrum.scanNumber, new List<Peptide>());
                    for (int i = currentSearchResultIdx; i < currentSearchResultIdx + TopN; i++)
                    {
                        result[spectrum.scanNumber].Add(Peptides[SearchResult[i]]);
                    }
                }
                else
                {
                    Console.WriteLine($"Warning: Found duplicate scan number {spectrum.scanNumber}. Skipping scan number...");
                }

                currentSearchResultIdx += TopN;
            }
        }

        /// <summary>
        /// Exports the results to csv format into a file with the given filename.
        /// </summary>
        /// <param name="filename">The output filename.</param>
        /// <returns>0 if the export was successful, 1 if was unsuccessful.</returns>
        public int export(string filename)
        {
            int status;

            try
            {
                var lines = new List<string>(){"ScanNumber;Peptides"};
                foreach (var item in result) {
                    var scanNr = item.Key;
                    var peptides = item.Value;
                    var line = scanNr.ToString() + ";";
                    foreach (var peptide in peptides)
                    {
                        line += (peptide.ToString() + ",");
                    }
                    lines.Add(line);
                }

                using (StreamWriter sw = new StreamWriter(filename))
                {
                    foreach (var line in lines)
                    {
                        sw.WriteLine(line);
                    }
                }

                status = 0;
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

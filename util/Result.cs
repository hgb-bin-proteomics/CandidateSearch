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
                        line += (peptide.toString() + ",");
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

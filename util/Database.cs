namespace CandidateSearch.util
{
    public class Peptide
    {
        public string sequence;
        public Dictionary<int, double> modifications;
        public bool isDecoy;

        public Peptide(string Sequence, Dictionary<int, double> Modifications, bool IsDecoy)
        {
            sequence = Sequence;
            modifications = Modifications;
            isDecoy = IsDecoy;
        }

        public int[] getEnconding(int massRange = 1300, int massMultiplier = 100)
        {
            var encoding = new List<int>();

            // code

            return encoding.ToArray();
        }

        public string toString()
        {
            string peptide = sequence + "[";
            foreach (var modification in modifications)
            {
                peptide = peptide + $"{modification.Key}:{modification.Value},";
            }

            if (isDecoy)
            {
                return "_" + peptide + "]";
            }

            return peptide + "]";
        }
    }

    public static class DatabaseReader
    {
        public static List<Peptide> readFASTA(string filename)
        {
            var peptides = new List<Peptide>();

            // code

            return peptides;
        }
    }
}

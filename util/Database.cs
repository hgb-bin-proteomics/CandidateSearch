using System.Text;

namespace CandidateSearch.util
{
    public class Peptide
    {
        public string sequence;
        public double mass;
        public Dictionary<int, double> modifications;
        public bool isDecoy;
        public List<double> ions;

        public Peptide(string Sequence, double Mass, Dictionary<int, double> Modifications, Settings IonSettings, bool IsDecoy)
        {
            sequence = Sequence;
            mass = Mass;
            modifications = Modifications;
            isDecoy = IsDecoy;
            ions = getIons(Sequence, Mass, IonSettings);
        }

        public int[] getEnconding(int massRange = 1300, int massMultiplier = 100)
        {
            var encoding = new List<int>();

            foreach (var ion in ions)
            {
                if (ion < massRange)
                {
                    encoding.Add((int) (ion * massMultiplier));
                }
            }

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

        private List<double> getIons(string Sequence, double Mass, Settings IonSettings)
        {
            var ions = new List<double>();

            for (int charge = 1; charge <= IonSettings.MAX_CHARGE; charge++)
            {
                double[] outIonsNoNL;
                MSAMANDA_IONCALCULATION.IonWithNL[] outIonsWithNL;
                MSAMANDA_IONCALCULATION.Modification[] mods = new MSAMANDA_IONCALCULATION.Modification[sequence.Length + 2];
                MSAMANDA_IONCALCULATION.IonCalculator.CalculateIons(out outIonsNoNL, 
                                                                    out outIonsWithNL, 
                                                                    Encoding.ASCII.GetBytes(Sequence),
                                                                    Mass,
                                                                    charge,
                                                                    mods,
                                                                    IonSettings.MAX_NEUTRAL_LOSSES,
                                                                    IonSettings.MAX_NEUTRAL_LOSS_MODS,
                                                                    0,
                                                                    1300,
                                                                    true,
                                                                    IonSettings.MAX_ALLOWED_CHARGE);

                ions.AddRange(outIonsNoNL);
            }

            return ions.Distinct().OrderBy(x => x).ToList();
        }
    }

    public static class DatabaseReader
    {
        public static List<Peptide> readFASTA(string filename, Settings settings, bool generateDecoys = false)
        {
            // digestion parameters set in method
            return MSAMANDA_FASTAPARSER.FASTAParser.DigestFasta(filename, settings, generateDecoys);
        }
    }
}

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

        public int[] getEnconding(int massRange = 5000, int massMultiplier = 100)
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

            double[] outIonsNoNL;
            MSAMANDA_IONCALCULATION.IonWithNL[] outIonsWithNL;
            MSAMANDA_IONCALCULATION.Modification[] mods = new MSAMANDA_IONCALCULATION.Modification[sequence.Length + 2];
            MSAMANDA_IONCALCULATION.IonCalculator.CalculateIons(out outIonsNoNL, 
                                                                out outIonsWithNL, 
                                                                sequence: Encoding.ASCII.GetBytes(Sequence),
                                                                mass: Mass,
                                                                charge: IonSettings.MAX_PRECURSOR_CHARGE,
                                                                mods: mods,
                                                                maxNumberNeutralLoss: IonSettings.MAX_NEUTRAL_LOSSES,
                                                                maxNumberNeutralLossModifications: IonSettings.MAX_NEUTRAL_LOSS_MODS,
                                                                lowerBound: 0,
                                                                upperBound: 5000,
                                                                mono: true,
                                                                maxAllowedChargeState: IonSettings.MAX_FRAGMENT_CHARGE);

            ions.AddRange(outIonsNoNL);

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

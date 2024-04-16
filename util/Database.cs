using System.Text;

namespace CandidateSearch.util
{
    /// <summary>
    /// Simplified peptide class that stores peptide/peptidoform information.
    /// </summary>
    public class Peptide
    {
        /// <summary>
        /// Amino acid sequence of the peptide.
        /// </summary>
        public string sequence { get; }
        /// <summary>
        /// Mass of the unmodified peptide.
        /// </summary>
        public double mass { get; }
        /// <summary>
        /// Dictionary mapping residue positions (0 based) to modification masses.
        /// </summary>
        public Dictionary<int, double> modifications { get; }
        /// <summary>
        /// Is the peptide a decoy peptide or target peptide.
        /// </summary>
        public bool isDecoy { get; }
        /// <summary>
        /// List of theoretical ion m/z values.
        /// </summary>
        public List<double> ions { get; }

        /// <summary>
        /// Constructor for a new peptide/pepidoform.
        /// </summary>
        /// <param name="Sequence">Sequence of amino acids.</param>
        /// <param name="Mass">Mass of the unmodified peptide.</param>
        /// <param name="Modifications">Dictionary that maps amino acid positions (0 based) to modification masses.</param>
        /// <param name="IonSettings">Settings used for ion calculation.</param>
        /// <param name="IsDecoy">Whether or not the peptide is a decoy peptide.</param>
        public Peptide(string Sequence, double Mass, Dictionary<int, double> Modifications, Settings IonSettings, bool IsDecoy)
        {
            sequence = Sequence;
            mass = Mass;
            modifications = Modifications;
            isDecoy = IsDecoy;
            ions = getIons(Sequence, Mass, Modifications, IonSettings);
        }

        /// <summary>
        /// Get the encoding vector of the peptide.
        /// </summary>
        /// <param name="massRange">Maximum m/z that should be considered while encoding. Has to match the specifications of VectorSearch.</param>
        /// <param name="massMultiplier">Precision of the encoding. Has to match the specifications of VectorSearch.</param>
        /// <returns>The encoding vector as an integer array.</returns>
        public int[] getEnconding(int massRange = 5000, int massMultiplier = 100)
        {
            var encoding = new List<int>();

            foreach (var ion in ions)
            {
                if (ion < massRange)
                {
                    encoding.Add((int) Math.Round(ion * massMultiplier));
                }
            }

            return encoding.Distinct().OrderBy(x => x).ToArray();
        }

        /// <summary>
        /// Constructs a string representation of the peptide.
        /// </summary>
        /// <returns>The string representation of the peptide.</returns>
        public override string ToString()
        {
            string peptide = sequence + "[";
            foreach (var modification in modifications)
            {
                peptide = peptide + $"{modification.Key}:{modification.Value}+";
            }

            peptide = peptide.TrimEnd(new char[] {'+'});

            if (isDecoy)
            {
                return "_" + peptide + "]";
            }

            return peptide + "]";
        }

        /// <summary>
        /// Adds a modification to the peptide if the peptide isn't already modified at that position.
        /// </summary>
        /// <param name="position">The position (0 based) of the modification.</param>
        /// <param name="mass">The modification mass.</param>
        /// <returns>True if the modification was added, false if there is already a modification on the specified residue.</returns>
        public bool addModification(int position, double mass)
        {
            if (modifications.ContainsKey(position))
                return false;

            modifications.Add(position, mass);
            return true;
        }

        /// <summary>
        /// Calculates theoretical ion m/z values for the peptide.
        /// </summary>
        /// <param name="Sequence">Sequence of amino acids.</param>
        /// <param name="Mass">Mass of the unmodified peptide.</param>
        /// <param name="Modifications">Dictionary that maps amino acid positions (0 based) to modification masses.</param>
        /// <param name="IonSettings">Settings used for ion calculation.</param>
        /// <returns>A list of theoretical ion m/z values.</returns>
        private List<double> getIons(string Sequence, double Mass, Dictionary<int, double> Modifications, Settings IonSettings)
        {
            var ions = new List<double>();

            double[] outIonsNoNL;
            MSAMANDA_IONCALCULATION.IonWithNL[] outIonsWithNL;
            MSAMANDA_IONCALCULATION.Modification[] mods = new MSAMANDA_IONCALCULATION.Modification[sequence.Length + 2];
            foreach (var mod in Modifications)
            {
                mods[mod.Key + 1] = new MSAMANDA_IONCALCULATION.Modification(title: mod.Key.ToString() + ":" + mod.Value.ToString(),
                                                                             name: mod.Key.ToString() + ":" + mod.Value.ToString(),
                                                                             mono: mod.Value,
                                                                             avg: mod.Value,
                                                                             aa: Sequence[mod.Key],
                                                                             fix: true, // this should technically be true/false depending on modification
                                                                             neutralLosses: new double[0],
                                                                             nTerminal: false,
                                                                             cTerminal: false,
                                                                             id: mod.Key.GetHashCode() + mod.Value.GetHashCode(),
                                                                             protein: false,
                                                                             maxOccurrence: 3);
            }
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

    /// <summary>
    /// Reader class to read a fasta file.
    /// </summary>
    public static class DatabaseReader
    {
        /// <summary>
        /// Reads and (tryptic) digests the given fasta file into a list of peptides.
        /// </summary>
        /// <param name="filename">The filename of the fasta file.</param>
        /// <param name="settings">Settings for digestion.</param>
        /// <param name="generateDecoys">Whether or not decoy peptides should be generated.</param>
        /// <returns>The list of peptides resulting from the digestion of the fasta file.</returns>
        public static List<Peptide> readFASTA(string filename, Settings settings, bool generateDecoys = false)
        {
            // digestion parameters set in method
            return MSAMANDA_FASTAPARSER.FASTAParser.DigestFasta(filename, settings, generateDecoys);
        }
    }
}

namespace MSAMANDA_CHEMICALUTILS
{
    public class ChemicalElement
    {
        public readonly string Name;
        public readonly string FullName;
        public readonly double MonoMass;
        public readonly double AvgMass;

        public ChemicalElement(string name, string fullName, double monoMass, double avgMass)
        {
            Name = name;
            FullName = fullName;
            MonoMass = monoMass;
            AvgMass = avgMass;
        }
    }

    public class AminoAcid
    {
        public readonly string Name;
        public readonly string ThreeLetterCode;
        public readonly char OneLetterCode;
        private readonly double AverageMass;
        private readonly double MonoisotopicMass;
        public readonly double ImmoniumMass;

        public double Mass(bool mono)
        {
            return mono ? MonoisotopicMass : AverageMass;
        }

        public AminoAcid(char oneLetter, string threeLetter, string name, double monoisotopicMass, double avgMass, double immoniumMass)
        {
            AverageMass = avgMass;
            MonoisotopicMass = monoisotopicMass;
            Name = name;
            ThreeLetterCode = threeLetter;
            OneLetterCode = oneLetter;
            ImmoniumMass = immoniumMass;
        }

        public AminoAcid(char oneLetter, string threeLetter, string name, double monoisotopicMass, double avgMass, bool mono = true)
        {
            AverageMass = avgMass;
            MonoisotopicMass = monoisotopicMass;
            Name = name;
            ThreeLetterCode = threeLetter;
            OneLetterCode = oneLetter;
            if (mono)
                ImmoniumMass = MonoisotopicMass + ChemicalUtils.Chemicals["H"].MonoMass -
                ChemicalUtils.Chemicals["C"].MonoMass - ChemicalUtils.Chemicals["O"].MonoMass;
            else
                ImmoniumMass = AverageMass + ChemicalUtils.Chemicals["H"].AvgMass -
                ChemicalUtils.Chemicals["C"].AvgMass - ChemicalUtils.Chemicals["O"].AvgMass;
        }
    }

    public class Permutations : IEquatable<Permutations>
    {
        public int N { get; set; }
        public int K { get; set; }
        public List<int[]> Elements { get; set; }

        public Permutations(int n, int k)
        {
            N = n;
            K = k;
            Elements = new List<int[]>();
            int[] data = new int[k];
            for (int i = 0; i < k; ++i)
                data[i] = i;
            Elements.Add(data);

            CalculatePermutations();
        }

        public Permutations(string key, List<int[]> elems)
        {
            string[] s = key.Split('_');
            N = Int32.Parse(s[0]);
            K = Int32.Parse(s[1]);
            Elements = elems;
        }

        public int[] Successor(int[] previous)
        {
            if (previous.Length == 0 || previous[0] == this.N - this.K)
                return null;

            int[] next = new int[K];

            int i;
            for (i = 0; i < this.K; ++i)
                next[i] = previous[i];

            for (i = this.K - 1; i > 0 && next[i] == this.N - this.K + i; --i) { }

            ++next[i];

            for (int j = i; j < this.K - 1; ++j)
                next[j + 1] = next[j] + 1;

            return next;
        }

        public void CalculatePermutations()
        {
            int[] last = Elements[0];
            int[] next = Successor(last);
            while (next != null)
            {
                Elements.Add(next);
                last = next;
                next = Successor(last);
            }
        }

        public bool Equals(Permutations other)
        {
            return (N == other.N && K == other.K);
        }
    }

    public static class ChemicalUtils
    {
        public static Dictionary<string, ChemicalElement> Chemicals = new Dictionary<string, ChemicalElement>() {
            {"H", new ChemicalElement("H", "Hydrogen", 1.007825035, 1.00794)},
            {"C", new ChemicalElement("C", "Carbon", 12.0, 12.0107)},
            {"N", new ChemicalElement("N", "Nitrogen", 14.003074, 14.0067)},
            {"O", new ChemicalElement("O", "Oxygen", 15.99491463, 15.9994)},
            {"p", new ChemicalElement("p", "Proton", 1.007276466812, 1.007276466812)},
            {"i", new ChemicalElement("i", "C13IsotopeDiff", 1.00335, 1.00335)}
        };

        public static Dictionary<char, AminoAcid> AminoAcids = new Dictionary<char, AminoAcid>() {
            {'A', new AminoAcid('A', "Ala", "Alanine", 71.037114, 71.0779)}, //71.03712
            {'R', new AminoAcid('R', "Arg", "Arginine", 156.101111, 156.1857)}, //156.10112
            {'N', new AminoAcid('N', "Asn", "Asparagine", 114.042927, 114.1026)}, //114.04293
            {'D', new AminoAcid('D', "Asp", "Aspartic acid", 115.026943, 115.0874)}, //115.02695
            {'C', new AminoAcid('C', "Cys", "Cysteine", 103.009185, 103.1429)}, //103.00919
            {'E', new AminoAcid('E', "Glu", "Glutamic acid", 129.042593, 129.114)}, //129.0426
            {'Q', new AminoAcid('Q', "Gln", "Glutamine", 128.058578, 128.1292)}, //128.05858
            {'G', new AminoAcid('G', "Gly", "Glycine", 57.021464, 57.0513)}, //57.02147
            {'H', new AminoAcid('H', "His", "Histidine", 137.058912, 137.1393)}, //137.05891
            {'I', new AminoAcid('I', "Ile", "Isoleucine", 113.084064, 113.1576)}, //113.08407
            {'L', new AminoAcid('L', "Leu", "Leucine", 113.084064, 113.1576)},
            {'K', new AminoAcid('K', "Lys", "Lysine", 128.094963, 128.1723)}, //128.09497
            {'M', new AminoAcid('M', "Met", "Methionine", 131.040485, 131.1961)}, //131.0405
            {'F', new AminoAcid('F', "Phe", "Phenylalanine", 147.068414, 147.1739)}, //147.06842
            {'P', new AminoAcid('P', "Pro", "Proline", 97.052764, 97.1152)}, //97.05277
            {'S', new AminoAcid('S', "Ser", "Serine", 87.032028, 87.0773)}, //87.03203
            {'T', new AminoAcid('T', "Thr", "Threonine", 101.047679, 101.1039)}, //101.04768
            {'U', new AminoAcid('U', "Sec", "Selenocysteine", 150.95363, 150.0379)},
            {'W', new AminoAcid('W', "Trp", "Tryptophane", 186.079313, 186.2099)}, //186.07932
            {'Y', new AminoAcid('Y', "Tyr", "Tyrosine", 163.06332, 163.1733)}, //163.06332
            {'V', new AminoAcid('V', "Val", "Valine", 99.068414, 99.1311)}, //99.06842
            {'J', new AminoAcid('J', "LeI", "Leu_Isoleu", 113.084064, 113.1576)},
            {'O', new AminoAcid('O', "Pyl", "Pyrrolysin", 237.14772677, 237.300713363271)},
            {'$', new AminoAcid('$', "", "N-Term", 0.0, 0.0)},
            {'^', new AminoAcid('^', "", "C-Term", 0.0, 0.0)}
        };


        public static Dictionary<string, Permutations> Combinations = new Dictionary<string, Permutations>();

        public static int GetMassIndex(double mass)
        {
            return (int)(mass * 10000);
        }

        public static double GetMassDoubleKey(double moverz)
        {
            return Math.Round(moverz, 7);
        }

        public static double CalculatePeptideMass(string sequence, bool mono)
        {
            double pepMass = Chemicals["H"].MonoMass * 2 + Chemicals["O"].MonoMass;
            if (!mono)
                pepMass = Chemicals["H"].AvgMass * 2 + Chemicals["O"].AvgMass;
            for (int i = 0; i < sequence.Length; ++i)
            {
                pepMass += AminoAcids[sequence[i]].Mass(mono);
            }
            return pepMass;
        }

        //public static double CalculatePeptideMass(string sequence, bool mono)
        //{
        //  return CalculatePeptideMass(sequence, 0, sequence.Length, mono);
        //}

        public static double CalculatePrecursorMassToCharge(int charge, double mass, bool mono)
        {
            if (mono)
                return (mass + charge * Chemicals["p"].MonoMass) / charge;
            else
                return (mass + charge * Chemicals["p"].AvgMass) / charge;
        }

        public static double CalculateUnchargedMass(double mass, int charge, bool mono)
        {
            if (mono)
                return ((mass * charge) - (charge) * Chemicals["p"].MonoMass);
            else
                return ((mass * charge) - (charge) * Chemicals["p"].AvgMass);
        }

        public static Permutations GetAscendingPermutations(int n, int k)
        {
            string key = n + "_" + k;
            lock (Combinations)
            {
                if (Combinations.ContainsKey(key))
                    return Combinations[key];
                else
                {
                    Permutations m = new Permutations(n, k);
                    lock (Combinations)
                    {
                        if (!Combinations.ContainsKey(key))
                        {
                            Combinations.Add(key, m);
                        }
                    }
                    return m;
                }
            }
        }

        public static double CalculateMassDiffInPPM(double expMass, double theoMass)
        {
            double ppm = (Math.Abs(expMass - theoMass) / theoMass) * (1000000.0);
            return ppm;
        }
    }
}

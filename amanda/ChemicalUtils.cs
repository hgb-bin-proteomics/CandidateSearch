/*
CODE FROM THE MS AMANDA SEARCH ENGINE
https://ms.imp.ac.at/?goto=msamanda

WHEN USED PLEASE CITE
https://doi.org/10.1021/pr500202e
https://doi.org/10.1002/rcm.9088

MS AMANDA LICENSE AGREEMENT

1.	This is the "MS Amanda" Freeware Software License Agreement, applying to the MS Amanda software with the copyright owners University of Applied Sciences Upper Austria GmbH (FH OÖ), Franz-Fritsch-Straße 11, 4600 Wels, and Research Institute of Molecular Pathology (IMP), Dr. Bohr-Gasse 7, 1030 Vienna, Austria. The software "MS Amanda" is provided as freeware. Authors of the software are Viktoria Dorfer, Marina Strobl, Sergey Maltsev, Peter Pichler, Stephan Winkler, and Karl Mechtler. All authors must be indicated on "MS Amanda". 
2.	Using the software MS Amanda you automatically agree to the terms and conditions contained within this Freeware Software License which is effective while using the software MS Amanda in the freeware version.
3.	The University of Applied Sciences Upper Austria GmbH (FH OÖ) and the Research Institute of Molecular Pathology (IMP) grant the Licensee a non-exclusive and non-transferable license to reproduce and use the software for personal purposes.
4.	Licensees may not rent, lease, lend, sell, reverse engineer, decompile, or disassemble "MS Amanda".
5.	Any proprietary rights, titles, ownership and intellectual property rights of "MS Amanda" remain with the University of Applied Sciences Upper Austria GmbH (FH OÖ) and the Research Institute of Molecular Pathology (IMP). The product is protected by copyright and intellectual property laws.
6.	The copyright holder reserves the right to reclassify this software as a non-freeware product at any time. This will not automatically affect the license agreement of previously distributed copies of the software. The University of Applied Sciences Upper Austria GmbH (FH OÖ) and the Research Institute of Molecular Pathology (IMP) may also terminate this agreement without any specified reason or if the licensee breaches any of its terms.
7.	Upon termination of this agreement the licensee is obliged to destroy all personal copies of "MS Amanda" obtained under this license.
8.	Licensees are not entitled to receive any documentation, technical support or updates to "MS Amanda".
9.	This agreement does not affect terms of other agreements regarding "MS Amanda".
10.	THE SOFTWARE "MS Amanda" IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
11.	Austrian Law is applicable to this license agreement. The relevant court for all disputes arising out of or in connection with this agreement is the competent court in Linz, Austria.


MS Amanda uses and incorporates third-party libraries or other resources 
that may be distributed under licenses different than the MS Amanda software.

MessagePack 
    https://msgpack.org
    Copyright (c) 2017 Yoshifumi Kawai and contributors
	MIT License

MIT License
=========================================

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
 */

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

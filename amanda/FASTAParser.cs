﻿/*
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

using System.Collections.Concurrent;
using System.Text.RegularExpressions;
using System.Text;
using MessagePack;
using CandidateSearch.util;
using MSAMANDA_CHEMICALUTILS;

namespace MSAMANDA_FASTAPARSER
{
    public class FASTAParser
    {
        public static List<Peptide> DigestFasta(string fastaFileName,
                                                Settings settings,
                                                bool generateDecoys = false, 
                                                double coreUsage = 0.75)
        {
            var trypsin = new Enzyme();
            trypsin.Name = "Trypsin";
            trypsin.CleavageSites = "KR";
            trypsin.CleavageInhibitors = "P";
            trypsin.Specificity = Enzyme.CLEAVAGE_SPECIFICITY.FULL;
            trypsin.Offset = 1;

            var proteins = ReadInFasta(fastaFileName, false);

            var peptides = DigestProteins(proteins, 
                                          trypsin, 
                                          settings,
                                          false, 
                                          coreUsage);

            if (generateDecoys)
            {
                var decoyProteins = new List<DBProtein>();
                foreach (var protein in proteins)
                {
                    var decoySequence = ReverseSequence(protein.Sequence);
                    var decoyProtein = new DBProtein(protein.DbProtRef.ProtIdentifier, -protein.DbProtRef.MappingId, decoySequence, true);
                    decoyProteins.Add(decoyProtein);
                }

                var decoyPeptides = DigestProteins(decoyProteins,
                                                   trypsin,
                                                   settings,
                                                   true,
                                                   coreUsage);

                peptides.AddRange(decoyPeptides);
            }

            return peptides;
        }

        private static string ReverseSequence(string seq)
        {
            char[] array = seq.ToCharArray();
            Array.Reverse(array);
            return new String(array);
        }

        private static List<Peptide> DigestProteins(List<DBProtein> proteins, 
                                                    Enzyme enzyme,
                                                    Settings settings,
                                                    bool isDecoy,
                                                    double coreUsage)
        {
            var opts = new ParallelOptions { MaxDegreeOfParallelism = (int) Math.Ceiling(Environment.ProcessorCount * coreUsage) };
            var concurrentPeptideList = new ConcurrentBag<List<DBPeptide>>();
            Parallel.ForEach(proteins, opts, (protein) => {
                var digester = new ProteinDigester(enzyme, settings.MAX_CLEAVAGES, true, settings.MIN_PEP_LENGTH, settings.MAX_PEP_LENGTH, protein);
                concurrentPeptideList.Add(digester.DigestProteinIntoList());
            });

            var peptideList = concurrentPeptideList.SelectMany(x => x).ToList();
            var massToPeptides = new DigesterDB();
            HelperMethods.MergeToDBDictionaries(peptideList, ref massToPeptides, opts);

            var peptides = new List<Peptide>();
            foreach (var item in massToPeptides.DbPeptidesDictMassKey)
            {
                var currentPeptides = item.Value;
                foreach (var peptide in currentPeptides)
                {
                    var peptidoforms = GetPeptidoforms(peptide, settings, isDecoy);
                    peptides.AddRange(peptidoforms);
                }
            }

            return peptides;
        }

        private static List<Peptide> GetPeptidoforms(DBPeptide dbPeptide, Settings settings, bool isDecoy)
        {
            var peptides = new List<Peptide>();
            var mods = new Dictionary<int, double>();

            if (settings.FIXED_MODIFICATIONS.Count > 0)
            {
                for (int i = 0; i < dbPeptide.Sequence.Length; i++)
                {
                    if (settings.FIXED_MODIFICATIONS.ContainsKey(dbPeptide.Sequence[i].ToString()))
                    {
                        mods.Add(i, settings.FIXED_MODIFICATIONS[dbPeptide.Sequence[i].ToString()]);
                    }
                }
            }

            var peptide = new Peptide(dbPeptide.Sequence, dbPeptide.Mass, mods, settings, isDecoy);
            peptides.Add(peptide);

            if (settings.VARIABLE_MODIFICATIONS.Count > 0)
            {
                foreach (var modification in settings.VARIABLE_MODIFICATIONS)
                {
                    addPeptidoformsForModification(peptides, modification, settings);
                }
            }

            return peptides;
        }

        private static void addPeptidoformsForModification(List<Peptide> peptides,
                                                           KeyValuePair<string, double> modification,
                                                           Settings settings)
        {
            var peptidoforms = new List<Peptide>();

            foreach (var peptide in peptides)
            {
                var possibleModificationSites = new List<int>();
                for (int i = 0; i < peptide.sequence.Length; i++)
                {
                    if (peptide.sequence[i].ToString() == modification.Key)
                    {
                        possibleModificationSites.Add(i);
                    }
                }

                var possibleCombinations = getAllPossibleCombinations(possibleModificationSites);

                foreach (var combination in possibleCombinations)
                {
                    var peptidoform = new Peptide(peptide.sequence,
                                                  peptide.mass,
                                                  new Dictionary<int, double>(),
                                                  settings,
                                                  peptide.isDecoy);

                    foreach (var mod in peptide.modifications)
                    {
                        peptidoform.addModification(mod.Key, mod.Value);
                    }


                    foreach (var position in combination)
                    {
                        peptidoform.addModification(position, modification.Value);
                    }

                    peptidoforms.Add(peptidoform);
                }
            }

            peptides.AddRange(peptidoforms);
        }

        private static List<List<int>> getAllPossibleCombinations(List<int> possibleModificationSites)
        {
            var possibleCombinations = new List<List<int>>();

            for (int i = 0; i < (1 << possibleModificationSites.Count); ++i)
            {
                var combination = new List<int>();
                for (int j = 0; j < possibleModificationSites.Count; ++j)
                {
                    if ((i & (1 << j)) != 0)
                    {
                        combination.Add(possibleModificationSites[j]);
                    }
                }
                possibleCombinations.Add(combination);
            }

            return possibleCombinations;
        }

        private static List<DBProtein> ReadInFasta(string fastaFileName, bool isDecoy)
        {
            List<DBProtein> proteins = new List<DBProtein>();
            string regexPatternSequence = "[^ARNDCEQGHILKMFPSTUWYVXBZJO]";
            int mappingID = 0;

            StreamReader fastaFileReader = new StreamReader(fastaFileName);
            try
            {
                string currentLine = "";
                string sequence = "";
                string identifier = "";

                while ((currentLine = fastaFileReader.ReadLine()) != null)
                {
                    if (currentLine.StartsWith(">", StringComparison.Ordinal))
                    {
                        if (!string.IsNullOrWhiteSpace(sequence))
                        {
                            //sequence = sequence.Replace('J', 'L');
                            if (Regex.IsMatch(sequence, regexPatternSequence))
                            {
                                var builder = new StringBuilder("Cannot parse ");
                                builder.Append("fasta file at identifier " + identifier + ". Sequence error. ");
                                Console.WriteLine(builder.ToString());
                                sequence = String.Empty;
                                identifier = string.Empty;
                                throw new Exception("Parsing error.");
                            }
                            proteins.Add(GenerateDbProtein(mappingID, isDecoy, identifier, sequence));
                            sequence = "";
                            mappingID++;
                        }

                        int index = currentLine.IndexOfAny(new[] { ' ', '|' }, 0);

                        if (index == -1)
                        {
                            identifier = currentLine.Substring(1);
                        }
                        else
                        {
                            identifier = currentLine.Substring(1, index - 1);
                            int indexNumber = identifier.IndexOfAny(new[]
                                {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'});
                            while (indexNumber == -1)
                            {
                                index = currentLine.IndexOfAny(new[] { ' ', '|' }, index + 1);
                                if (index == -1)
                                {
                                    identifier = currentLine.Substring(1);
                                    break;
                                }

                                identifier = currentLine.Substring(1, index - 1);
                                indexNumber = identifier.IndexOfAny(new[]
                                    {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'});
                            }
                        }
                    }
                    else
                    {
                        sequence += currentLine.ToUpper();
                    }
                }

                if (!string.IsNullOrWhiteSpace(sequence))
                {
                    if (Regex.IsMatch(sequence, regexPatternSequence))
                    {
                        var builder = new StringBuilder("Cannot parse ");
                        builder.Append("fasta file at identifier " + identifier + ". Sequence error. ");
                        Console.WriteLine(builder.ToString());
                        throw new Exception("Parsing error.");
                    }
                    else
                    {
                        proteins.Add(GenerateDbProtein(mappingID, isDecoy, identifier, sequence));
                        sequence = "";
                        mappingID++;
                    }
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("Fasta file error");
                Console.WriteLine(e.ToString());
                throw;
            }
            finally
            {
                fastaFileReader.Close();
            }

            return proteins;
        }

        private static DBProtein GenerateDbProtein(int mappingID, bool isDecoy, string identifier, string sequence)
        {

            if (isDecoy)
            {
                //mark proteins as decoys by negative ProteinIDs;
                return new DBProtein(identifier, -mappingID, sequence);
            }
            return new DBProtein(identifier, mappingID, sequence);
        }
    }

    public class DBProtein
    {
        public string Sequence { get; set; }
        public bool IsDecoy { get; set; }
        public DBProtRef DbProtRef { get; set; }

        public DBProtein(string identifier, int id, string sequence, bool isDecoy = false)
        {
            var ok = int.TryParse(identifier, out int identy);
            if (ok)
            {
                DbProtRef = new DBProtRef
                {
                    ProtId = identy,
                    MappingId = id,
                    ProtIdentifier = identifier
                };
            }
            else
            {
                DbProtRef = new DBProtRef
                {
                    ProtId = id,
                    MappingId = id,
                    ProtIdentifier = identifier
                };
            }
            Sequence = sequence;
            IsDecoy = isDecoy;
        }
    }

    [MessagePackObject]
    public struct DBProtRef
    {
        [Key(0)]
        public int ProtId { get; set; }
        [Key(1)]
        public int MappingId { get; set; }
        [Key(2)]
        public string ProtIdentifier { get; set; }

    }

    public class Enzyme
    {
        public Enzyme()
        {
            CleavageInhibitors = "";
            CleavageSites = "";
            Offset = 0;
            Specificity = CLEAVAGE_SPECIFICITY.FULL;
            Name = "";
        }

        public enum CLEAVAGE_SPECIFICITY { FULL, SEMI, SEMI_N, SEMI_C };

        public string CleavageSites { get; set; }
        public string CleavageInhibitors { get; set; }
        public CLEAVAGE_SPECIFICITY Specificity { get; set; }
        public int Offset { get; set; }
        public string Name { get; set; }

        public Regex TheRegex => regexInit.Value;
        private Lazy<Regex> regexInit => new Lazy<Regex>(() =>
        {
            if (Offset == 1)
            {
                //C-Terminal
                StringBuilder builder = new StringBuilder("(?<=[");
                builder.Append(CleavageSites);
                builder.Append("])");
                if (!String.IsNullOrEmpty(CleavageInhibitors))
                {
                    builder.Append("(?=[^");
                    builder.Append(CleavageInhibitors);
                    builder.Append("])");
                }
                // "(?<=[" + Enzyme.CleavageSites + "])(?=[^" + Enzyme.CleavageInhibitors + "])"
                return new Regex(builder.ToString());
            }
            else
            {
                //N-Terminal
                StringBuilder builder = new StringBuilder("(?=[");
                builder.Append(CleavageSites);
                builder.Append("]");
                if (!String.IsNullOrEmpty(CleavageInhibitors))
                {
                    builder.Append("[^");
                    builder.Append(CleavageInhibitors);
                    builder.Append("]");
                }
                builder.Append(")");
                // "(?=[" + Enzyme.CleavageSites + "][^" + Enzyme.CleavageInhibitors + "])"
                return new Regex(builder.ToString());
            }
        });
    }

    [MessagePackObject]
    public class DBPeptide
    {
        [Key(0)]
        public double Mass { get; set; }
        [Key(1)]
        public int MassInt { get; set; }
        [Key(2)]
        public string Sequence { get; set; }
        [Key(3)]
        public bool ProteinStartFlag { get; set; }
        [Key(4)]
        public List<DBProtRef> DbProtRefs { get; set; }
        [Key(5)]
        public string SequenceOriginal { get; set; }
        [Key(6)]
        public int MissedCleavages { get; set; }
        [IgnoreMember]
        public int SeqHash { get; set; }

        public DBPeptide()
        { 
            Sequence = ""; 
            DbProtRefs = new List<DBProtRef>(); 
            SequenceOriginal = "";
        }

        public DBPeptide(string sequence, string sequenceOriginal, int missedCleavages, bool proteinStartFlag, DBProtRef protRef)
        {
            MissedCleavages = missedCleavages;
            Sequence = sequence;
            SequenceOriginal = sequenceOriginal;
            ProteinStartFlag = proteinStartFlag;
            SeqHash = CreateMD5();

            DbProtRefs = new List<DBProtRef>(1) { protRef };
        }

        internal void AddToProtRefs(DBPeptide pep)
        {
            if (DbProtRefs != null && DbProtRefs.Count > 0)
            {
                DbProtRefs = DbProtRefs.Concat(pep.DbProtRefs).Distinct().ToList();
            }
            else
            {
                DbProtRefs = pep.DbProtRefs;
            }
        }

        public int CreateMD5()
        {
            return (Sequence, MissedCleavages, ProteinStartFlag).GetHashCode();
        }
    }

    public static class HelperMethods
    {
        public static DBPeptide FindPeptideWithSameHash(this List<DBPeptide> peptides, DBPeptide pep)
        {
            for (var i = 0; i < peptides.Count; i++)
            {
                var x = peptides[i];
                if (x.SeqHash == pep.SeqHash)
                {
                    return x;
                }
            }
            return null;
        }

        public static bool IsBetweenExcludeBounds(this double target, double start, double end)
        {
            return target > start && target < end;
        }

        public static bool IsBetweenIncludeBounds(this int target, int start, int end)
        {
            return target >= start && target <= end;
        }

        public static void MergeToDBDictionaries(List<DBPeptide> dbFrom, ref DigesterDB dbTo, ParallelOptions parallelOptions)
        {
            dbTo.DbPeptidesDictMassKey = dbFrom.GroupBy(x => x.MassInt)
                                               .AsParallel()
                                               .WithDegreeOfParallelism(parallelOptions.MaxDegreeOfParallelism)
                                               .Select(g => {
                var theList = new List<DBPeptide>();
                foreach (var pep in g)
                {
                    var itemWithSameHash = theList.FindPeptideEqualTo(pep);
                    //var itemWithSameHash = theList.Find(x => DBPeptideEquals(x, pep));
                    if (itemWithSameHash != null)
                    {
                        itemWithSameHash.AddToProtRefs(pep);
                    }
                    else
                    {
                        theList.Add(pep);
                    }
                }
                return (g.Key, theList);
            }).ToDictionary(t => t.Key, t => t.theList);
        }

        private static DBPeptide FindPeptideEqualTo(this List<DBPeptide> list, DBPeptide comparePep)
        {
            for (var i = 0; i < list.Count; i++)
            {
                var x = list[i];
                if (DBPeptideEquals(comparePep, x))
                {
                    return x;
                }
            }
            return null;
        }

        private static bool DBPeptideEquals(DBPeptide x, DBPeptide pep)
        {
            if (x.SeqHash != pep.SeqHash)
            {
                return false;
            }
            if (x.ProteinStartFlag != pep.ProteinStartFlag)
            {
                return false;
            }
            if (x.MissedCleavages != pep.MissedCleavages)
            {
                return false;
            }
            if (x.Sequence != pep.Sequence)
            {
                return false;
            }
            return true;
        }
    }

    public class ProteinDigester
    {
        private readonly Enzyme _enzyme;
        private readonly int _missedCleavages;
        private readonly bool _useMonoisotopicMass;
        private readonly int _minPepLength;
        private readonly int _maxPepLength;
        private readonly DBProtein _dbProtein;

        private List<DBPeptide> _dbPeptides;

        struct PeptideInfo
        {
            public int Start;
            public int MissedCleavages;
        };

        public ProteinDigester(Enzyme enzyme, int missedCleavages, bool useMonoisotopicMass, int minPepLength, int maxPepLength, DBProtein dbProtein)
        {
            _enzyme = enzyme;
            _missedCleavages = missedCleavages;
            _useMonoisotopicMass = useMonoisotopicMass;
            _minPepLength = minPepLength;
            _maxPepLength = maxPepLength;
            _dbProtein = dbProtein;
            _dbPeptides = new List<DBPeptide>();
        }

        public List<DBPeptide> DigestProteinIntoList()
        {
            _dbPeptides = new List<DBPeptide>();
            DigestSingleProtein(_dbProtein);
            return _dbPeptides;
        }

        private void DigestSingleProtein(DBProtein protein)
        {

            if (_enzyme.CleavageSites == "X")
            {
                //no-enzyme search
                DigestSingleProteinWithNoEnzyme(protein);
            }
            else if (_enzyme.CleavageSites == "")
            {
                //no cleavage (peptide db or top down)
                SaveSinglePeptide(new DBPeptide(protein.Sequence, protein.Sequence, 0, IsProteinStart(0, protein.Sequence[0]), protein.DbProtRef));
            }
            else
            {
                var peptides = SplitProtein(_enzyme, protein, _missedCleavages, 0);
                foreach (var pep in peptides)
                {
                    SaveSinglePeptide(pep);
                }

                if (protein.Sequence.StartsWith("M", StringComparison.Ordinal) && _enzyme.CleavageSites != "X" && _enzyme.CleavageSites != "")
                {
                    var prot = new DBProtein(protein.DbProtRef.ProtIdentifier, protein.DbProtRef.MappingId, protein.Sequence.Substring(1));
                    var otherPeptides = SplitProtein(_enzyme, prot, _missedCleavages, 1);

                    foreach (var pep in otherPeptides)
                    {
                        var peptideWithSameHash = peptides.FindPeptideWithSameHash(pep);
                        if (peptideWithSameHash == null)
                        {
                            SaveSinglePeptide(pep);
                            continue;
                        }
                        if (peptideWithSameHash.MissedCleavages == pep.MissedCleavages) continue;

                        SaveSinglePeptide(pep);
                    }
                }
            }
        }

        private void DigestSingleProteinWithNoEnzyme(DBProtein protein)
        {
            SortedSet<int> positionsOfX = new SortedSet<int>();
            SortedSet<int> positionsOfZ = new SortedSet<int>();
            SortedSet<int> positionsOfB = new SortedSet<int>();
            MatchCollection collX = Regex.Matches(protein.Sequence, "X");
            MatchCollection collZ = Regex.Matches(protein.Sequence, "Z");
            MatchCollection collB = Regex.Matches(protein.Sequence, "B");

            for (int m = 0; m < collX.Count; ++m)
            {
                positionsOfX.Add(collX[m].Index);
            }

            for (int m = 0; m < collZ.Count; ++m)
            {
                positionsOfZ.Add(collZ[m].Index);
            }

            for (int m = 0; m < collB.Count; ++m)
            {
                positionsOfB.Add(collB[m].Index);
            }

            for (int s = 0; s < protein.Sequence.Length; ++s)
            {
                for (int l = 1; l < protein.Sequence.Length - s + 1; ++l)
                {
                    var seqL = l - s;
                    if (seqL < _minPepLength)
                        continue;
                    if (seqL > _maxPepLength)
                        break;
                    if (l > 150)
                        break;
                    SortedSet<int> xInSeq = positionsOfX.GetViewBetween(s, s + l);
                    if (xInSeq.Count > 1)
                        break;
                    SortedSet<int> bInSeq = positionsOfB.GetViewBetween(s, s + l);
                    if (bInSeq.Count > 3)
                        break;
                    SortedSet<int> zInSeq = positionsOfZ.GetViewBetween(s, s + l);
                    if (zInSeq.Count > 3 || zInSeq.Count + bInSeq.Count > 3)
                        break;
                    if (xInSeq.Count > 0 || bInSeq.Count > 0 || zInSeq.Count > 0)
                    {
                        //handle replacement characters
                        SaveSinglePeptide(new DBPeptide(protein.Sequence.Substring(s, l), protein.Sequence.Substring(s, l), 0,
                          IsProteinStart(s, protein.Sequence[0]), protein.DbProtRef));
                    }
                    else
                    {
                        //save pep
                        string sseq = protein.Sequence.Substring(s, l);
                        double mass = ChemicalUtils.CalculatePeptideMass(sseq, _useMonoisotopicMass);

                        if (mass.IsBetweenExcludeBounds(200, 6000))
                        //if (!(mass > 200) || !(mass < 6000))
                        //{
                        //  if (mass > 6000)
                        //    break;
                        //}
                        //else
                        {
                            SavePeptide(new DBPeptide(sseq, sseq, 0, IsProteinStart(s, protein.Sequence[0]), protein.DbProtRef));
                        }
                    }
                }
            }
        }

        private void SaveSinglePeptide(DBPeptide pep)
        {
            if (String.IsNullOrEmpty(pep.Sequence))
                return;
            if (pep.Sequence.Length < _minPepLength || pep.Sequence.Length > _maxPepLength)
                return;

            if (pep.Sequence.Contains('B') || pep.Sequence.Contains('Z') || pep.Sequence.Contains('X'))
            {
                string temp = pep.Sequence.Replace("X", "");
                if (pep.Sequence.Length - temp.Length > 1)
                {
                    return;
                }

                temp = pep.Sequence.Replace("B", "");
                temp = temp.Replace("Z", "");
                if (pep.Sequence.Length - temp.Length > 3)
                {
                    return;
                }

                List<char> replacements = new List<char>();
                for (int i = 0; i < pep.Sequence.Length; ++i)
                {
                    replacements.Add('#');
                }

                if (pep.Sequence.Contains('X'))
                {
                    GenerateCombinationsForX(pep, replacements);
                }
                else
                {
                    CalculateMass(pep, replacements);
                }
            }
            else
            {
                pep.Mass = ChemicalUtils.CalculatePeptideMass(pep.Sequence, _useMonoisotopicMass);
                SavePeptide(pep);
            }
        }

        private static List<DBPeptide> SplitProtein(Enzyme enzyme, DBProtein prot, int missedCleavages, int offset)
        {
            string[] peptides = enzyme.TheRegex.Split(prot.Sequence);
            List<DBPeptide> peps = new List<DBPeptide>(offset == 0 ? peptides.Length * 3 : 3); // guessing magic numbers
            int start = offset;

            Dictionary<string, PeptideInfo> allPeps = new Dictionary<string, PeptideInfo>();
            for (int i = 0; i < peptides.Length; ++i)
            {
                if (offset > 0 && i > 0)
                    break;
                if (!String.IsNullOrEmpty(peptides[i]))
                {
                    string myPeptide = peptides[i];
                    if (!allPeps.ContainsKey(myPeptide))
                    {
                        //peps.Add(new DBPeptide(myPeptide, myPeptide, 0, IsProteinStart(0, myPeptide[0]), prot.DbProtRef));
                        peps.Add(new DBPeptide(myPeptide, myPeptide, 0, IsProteinStart(myPeptide, prot.Sequence), prot.DbProtRef));
                        allPeps.Add(peptides[i], new PeptideInfo
                        {
                            Start = start,
                            MissedCleavages = 0
                        });
                    }

                    int j = 1;
                    string pep = peptides[i];

                    int count = 0;
                    while (count < missedCleavages)
                    {
                        if (i + j < peptides.Length)
                        {
                            if (!String.IsNullOrEmpty(peptides[i + j]))
                            {
                                pep += peptides[i + j];
                                if (!allPeps.ContainsKey(pep))
                                {
                                    //peps.Add(new DBPeptide(pep, pep, count + 1, IsProteinStart(0, pep[0]), prot.DbProtRef));
                                    peps.Add(new DBPeptide(pep, pep, count + 1, IsProteinStart(pep, prot.Sequence), prot.DbProtRef));
                                    allPeps.Add(pep, new PeptideInfo
                                    {
                                        Start = start,
                                        MissedCleavages = count + 1
                                    });
                                }

                                ++count;
                            }

                            ++j;
                        }
                        else
                        {
                            break;
                        }
                    }

                    start += peptides[i].Length;
                }
            }

            if (enzyme.CleavageSites != "X" && (enzyme.Specificity == Enzyme.CLEAVAGE_SPECIFICITY.SEMI ||
                                                enzyme.Specificity == Enzyme.CLEAVAGE_SPECIFICITY.SEMI_C ||
                                                enzyme.Specificity == Enzyme.CLEAVAGE_SPECIFICITY.SEMI_N))
            {
                foreach (string peptide in allPeps.Keys.ToArray())
                {
                    PeptideInfo info = allPeps[peptide];
                    for (int i = 1; i < peptide.Length - 1; ++i)
                    {
                        string newFront = peptide.Substring(i);
                        string newBack = peptide.Substring(0, peptide.Length - i);
                        if (!allPeps.ContainsKey(newFront) && (enzyme.Specificity == Enzyme.CLEAVAGE_SPECIFICITY.SEMI ||
                                                               enzyme.Specificity == Enzyme.CLEAVAGE_SPECIFICITY.SEMI_C))
                        {
                            //peps.Add(new DBPeptide(newFront, newFront, info.MissedCleavages, IsProteinStart(0, newFront[0]), prot.DbProtRef));
                            peps.Add(new DBPeptide(newFront, newFront, info.MissedCleavages, IsProteinStart(newFront, prot.Sequence), prot.DbProtRef));

                            allPeps.Add(newFront, new PeptideInfo
                            {
                                Start = info.Start + i,
                                MissedCleavages = info.MissedCleavages
                            });
                        }

                        if (allPeps.ContainsKey(newBack) || (enzyme.Specificity != Enzyme.CLEAVAGE_SPECIFICITY.SEMI &&
                                                             enzyme.Specificity != Enzyme.CLEAVAGE_SPECIFICITY.SEMI_N)) continue;

                        //peps.Add(new DBPeptide(newBack, newBack, info.MissedCleavages, IsProteinStart(0, newBack[0]), prot.DbProtRef));
                        peps.Add(new DBPeptide(newBack, newBack, info.MissedCleavages, IsProteinStart(newBack, prot.Sequence), prot.DbProtRef));
                        allPeps.Add(newBack, new PeptideInfo
                        {
                            Start = info.Start,
                            MissedCleavages = info.MissedCleavages
                        });
                    }
                }
            }

            return peps;
        }

        private void SavePeptide(DBPeptide pep)
        {
            if (pep.Sequence.Length < 2)
                return;
            if (!pep.Sequence.Length.IsBetweenIncludeBounds(_minPepLength, _maxPepLength))
                return;

            pep.Mass = ChemicalUtils.CalculatePeptideMass(pep.Sequence, _useMonoisotopicMass);
            pep.MassInt = (int)(pep.Mass);
            pep.SeqHash = pep.CreateMD5();
            _dbPeptides.Add(pep);
        }

        private DBPeptide CreatePeptide(string sequence, string origSequence, int missedCleavages, bool proteinStartFlag, DBProtRef protRef)
        {
            var pep = new DBPeptide(sequence, origSequence, missedCleavages, proteinStartFlag, protRef);
            if (pep.Sequence.Length < 2)
                return null;
            if (pep.Sequence.Length < _minPepLength || pep.Sequence.Length > _maxPepLength)
                return null;
            pep.Mass = ChemicalUtils.CalculatePeptideMass(pep.Sequence, _useMonoisotopicMass);
            pep.MassInt = (int)(pep.Mass);
            pep.SeqHash = pep.CreateMD5();
            return pep;
        }

        private static bool IsProteinStart(int startPosition, char protStarter)
        {
            return (startPosition == 0 || (startPosition == 1 && protStarter == 'M'));
        }

        private static bool IsProteinStart(string pep, string prot)
        {
            return IsProteinStart(prot.IndexOf(pep), prot[0]);
        }

        private void CalculateMass(DBPeptide info,
          List<char> replacements)
        {
            if (info.Sequence.Contains('B'))
            {
                int firstIndex = info.Sequence.IndexOf('B');
                string changedPeptide = info.Sequence;
                char[] temp = new char[replacements.Count];
                replacements.CopyTo(temp);
                List<char> currentReplacements = new List<char>(temp)
                {
                    [firstIndex] = 'N'
                };
                changedPeptide = changedPeptide.Insert(firstIndex, "N");
                changedPeptide = changedPeptide.Remove(firstIndex + 1, 1);
                info.Sequence = changedPeptide;
                CalculateMass(info, currentReplacements);
                changedPeptide = changedPeptide.Insert(firstIndex, "D");
                changedPeptide = changedPeptide.Remove(firstIndex + 1, 1);
                currentReplacements[firstIndex] = 'D';
                info.Sequence = changedPeptide;
                CalculateMass(info, currentReplacements);
            }
            else if (info.Sequence.Contains('Z'))
            {
                int firstIndex = info.Sequence.IndexOf('Z');
                string changedPeptide = info.Sequence;
                char[] temp = new char[replacements.Count];
                replacements.CopyTo(temp);
                List<char> currentReplacements = new List<char>(temp)
                {
                    [firstIndex] = 'E'
                };
                changedPeptide = changedPeptide.Insert(firstIndex, "E");
                changedPeptide = changedPeptide.Remove(firstIndex + 1, 1);
                info.Sequence = changedPeptide;
                CalculateMass(info, currentReplacements);
                changedPeptide = changedPeptide.Insert(firstIndex, "Q");
                changedPeptide = changedPeptide.Remove(firstIndex + 1, 1);
                currentReplacements[firstIndex] = 'Q';
                info.Sequence = changedPeptide;
                CalculateMass(info, currentReplacements);
            }
            else
            {
                SavePeptide(info);
                //FinalCalculation(info);
            }
        }

        private void FinalCalculation(DBPeptide info)
        {
            double mass = ChemicalUtils.CalculatePeptideMass(info.Sequence, _useMonoisotopicMass);
            SavePeptide(info);
        }

        private void GenerateCombinationsForX(DBPeptide info,
          List<char> replacements)
        {
            foreach (char aa in ChemicalUtils.AminoAcids.Keys)
            {
                if (aa != '^' && aa != '$' && aa != 'J')
                {
                    int firstIndex = info.Sequence.IndexOf('X');
                    string changedPeptide = info.Sequence;
                    changedPeptide = changedPeptide.Insert(firstIndex, aa.ToString());
                    changedPeptide = changedPeptide.Remove(firstIndex + 1, 1);
                    //char[] temp = new char[replacements.Count];
                    //replacements.CopyTo(temp);

                    List<char> currentReplacement = new List<char>(replacements.ToArray())
                    {
                        [firstIndex] = aa
                    };
                    //string currentReplacement = variableAAs + aa.ToString();

                    if (changedPeptide.Contains('X'))
                    {
                        //should not occur
                    }
                    else
                    {
                        var newPep = new DBPeptide
                        {
                            MissedCleavages = info.MissedCleavages,
                            Sequence = changedPeptide,
                            SequenceOriginal = info.SequenceOriginal,
                            ProteinStartFlag = info.ProteinStartFlag,
                            DbProtRefs = info.DbProtRefs
                        };
                        CalculateMass(newPep, currentReplacement);
                    }
                }
            }
        }
    }

    [MessagePackObject]
    public class DigesterDB
    {
        [Key(0)]
        public Dictionary<int, List<DBPeptide>> DbPeptidesDictMassKey { get; set; }

        public DigesterDB()
        {
            DbPeptidesDictMassKey = new Dictionary<int, List<DBPeptide>>();
        }
    }
}

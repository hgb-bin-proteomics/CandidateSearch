using System.Text.RegularExpressions;
using System.Text;
using MessagePack;

namespace MSAMANDA_FASTAPARSER
{
    public class FASTAParser
    {
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
}

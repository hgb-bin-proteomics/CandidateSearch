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

using System.Globalization;
using System.Text.RegularExpressions;

namespace MSAMANDA_MGFPARSER
{
    public class AMassCentroid
    {
        public double Position { get; set; }
        public double Intensity { get; set; }
        public int Charge { get; set; }

        public AMassCentroid Clone()
        {
            return new AMassCentroid() { Position = this.Position, Intensity = this.Intensity, Charge = this.Charge };
        }

        public AMassCentroid Clone(int charge)
        {
            return new AMassCentroid() { Position = this.Position, Intensity = this.Intensity, Charge = charge };
        }
    }

    public class Precursor
    {
        public double MOverZ { get; private set; }
        public double UnChargedMass { get; private set; }
        public double Intensity { get; set; }
        public int Charge { get; private set; }
        public int Rank { get; set; }

        public void SetMassCharge(double mz, int charge, bool mono)
        {
            MOverZ = mz;
            Charge = charge;
            UnChargedMass = MSAMANDA_CHEMICALUTILS.ChemicalUtils.CalculateUnchargedMass(MOverZ, Charge, mono);
        }
    }

    public class Spectrum
    {
        public int ScanNumber { get; set; }
        public int SpectrumId { get; set; }
        public double RT { get; set; }
        public List<AMassCentroid> FragmentsPeaks { get; set; }
        public Dictionary<int, double> ImmunePeaks { get; set; }
        public SortedSet<double> ImmuneMasses { get; set; }
        public Precursor Precursor { get; set; }

        public Spectrum()
        {
            ScanNumber = 0;
            RT = 0.0;
            SpectrumId = 0;
            Precursor = new Precursor();
        }
    }

    public static class MGFParser
    {
        public static List<Spectrum> ParseNextSpectra(string filename)
        {
            List<Spectrum> spectra = new List<Spectrum>();
            int nrOfReadSpectra = 0;
            int nrOfSpectra = 0;
            double mOverZ = 0.0;
            int charge = 0;
            Dictionary<int, double> peaks = new Dictionary<int, double>();
            Spectrum currentSpectrum = null;
            bool isCorrect = true;
            string lastScannumber = string.Empty;
            string title = "";
            int _scanId = 0;
            var SpectTitleMap = new Dictionary<int, string>();

            try
            {
                using (StreamReader sr = new StreamReader(filename))
                {
                    while (!sr.EndOfStream)
                    {
                        var line = sr.ReadLine();

                        if (line == null) continue;

                        // error in parsing spectrum, search for next one
                        if (!isCorrect && (!line.ToUpper().StartsWith("BEGIN IONS", StringComparison.Ordinal)))
                        {
                            continue;
                        }

                        if (line.ToUpper().StartsWith("BEGIN IONS", StringComparison.Ordinal))
                        {
                            if (currentSpectrum != null && isCorrect)
                            {
                                Console.WriteLine(
                                    "  Skipping spectrum with scannumber '" + currentSpectrum.ScanNumber +
                                    "'. No end ions found."
                                    );
                            }

                            ++nrOfSpectra;
                            isCorrect = true;
                            peaks = new Dictionary<int, double>();
                            mOverZ = 0.0;
                            charge = 0;
                            currentSpectrum = new Spectrum
                            {
                                FragmentsPeaks = new List<AMassCentroid>(),
                                ImmuneMasses = new SortedSet<double>(),
                                ImmunePeaks = new Dictionary<int, double>(),
                                SpectrumId = _scanId
                            };
                        }
                        else if (LineCanBeIgnored(line)) { }
                        else
                        {
                            if (currentSpectrum == null)
                            {
                                isCorrect = false;
                                Console.WriteLine(
                                    "  Skipping spectrum after scannumber '" + lastScannumber + "'. No begin ions found."
                                    );
                                continue;
                            }

                            if (line.ToUpper().StartsWith("TITLE", StringComparison.Ordinal))
                            {
                                int inx = line.IndexOf("=");
                                title = line.Substring(inx + 1);
                                if (line.ToUpper().Contains("SCAN", StringComparison.Ordinal))
                                {
                                    Match s = Regex.Match(line.ToUpper(), @"SCAN.?[:=\s]\s?([0-9]+)");
                                    if (s.Success)
                                        currentSpectrum.ScanNumber = Int32.Parse(s.Groups[1].Value);
                                }
                                else if (line.ToUpper().Contains("INDEX", StringComparison.Ordinal))
                                {
                                    Match s = Regex.Match(line.ToUpper(), @"INDEX.?[:=\s]\s?([0-9]+)");
                                    if (s.Success)
                                        currentSpectrum.ScanNumber = Int32.Parse(s.Groups[1].Value);
                                }
                            }
                            else if (line.ToUpper().StartsWith("PEPMASS", StringComparison.Ordinal))
                            {
                                mOverZ = ParseMOverZ(line);
                            }
                            else if (line.ToUpper().StartsWith("CHARGE", StringComparison.Ordinal))
                            {
                                charge = ParseCharge(line);
                            }
                            else if (line.ToUpper().StartsWith("RTINSECONDS", StringComparison.Ordinal))
                            {
                                currentSpectrum.RT = ParseRt(line);
                            }
                            else if (line.ToUpper().StartsWith("SCANS", StringComparison.Ordinal))
                            {
                                currentSpectrum.ScanNumber = ParseScanNumber(line);
                            }
                            else if (line.ToUpper().StartsWith("END IONS", StringComparison.Ordinal))
                            {
                                if (currentSpectrum.ScanNumber == 0)
                                {
                                    if (string.IsNullOrEmpty(title))
                                    {
                                        Console.WriteLine(
                                            "  Skipping spectrum after scannumber '" +
                                            lastScannumber + "'. No title or scan number found."
                                            );
                                        isCorrect = false;
                                        continue;
                                    }

                                    string[] titleArr = title.Split('.');
                                    if (titleArr.Length > 3)
                                        currentSpectrum.ScanNumber =
                                            Int32.Parse(titleArr[titleArr.Length - 3]); // vorvorletztes item; => firstScan
                                }

                                if (currentSpectrum.FragmentsPeaks.Count == 0 || mOverZ == 0)
                                {
                                    string text = title;
                                    if (currentSpectrum.ScanNumber != 0)
                                        text = currentSpectrum.ScanNumber.ToString();

                                    if (mOverZ == 0)
                                    {
                                        Console.WriteLine(
                                            "  Skipping spectrum with scannumber '" + text + "'. No mass value found."
                                            );
                                    }
                                    else
                                    {
                                        Console.WriteLine(
                                            "  Skipping spectrum with scannumber '" + text + "'. No peaks found."
                                            );
                                    }

                                    isCorrect = false;
                                    continue;
                                }

                                // masses and peaks not needed, new calculation from fragmentpeaks
                                // masses and peaks are overwritten in Spectrum.PrepareForSearch()

                                if (charge == 0)
                                {
                                    foreach (int consideredCharge in new List<int>() { 2, 3, 4, 5, 6 })
                                    {
                                        var s = GenerateSpectrum(currentSpectrum.FragmentsPeaks,
                                            currentSpectrum.ScanNumber, _scanId, currentSpectrum.RT, mOverZ,
                                            consideredCharge);
                                        spectra.Add(s);
                                        _scanId++;
                                    }
                                }
                                else
                                {
                                    currentSpectrum.Precursor.SetMassCharge(mOverZ, charge, true);
                                    spectra.Add(currentSpectrum);
                                    SpectTitleMap.Add(_scanId, title.Trim());
                                    _scanId++;
                                }


                                ++nrOfReadSpectra;
                                lastScannumber = currentSpectrum.ScanNumber.ToString();
                                title = string.Empty;
                                peaks = new Dictionary<int, double>();
                                mOverZ = 0.0;
                                charge = 0;
                                currentSpectrum = null;
                            }
                            // fragments
                            else
                            {
                                string[] parts = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                                if (parts.Length == 1)
                                {
                                    parts = line.Split(new[] { '\t' }, StringSplitOptions.RemoveEmptyEntries);
                                }

                                if (parts.Length == 2 || parts.Length == 3)
                                {
                                    double mass;
                                    double intensity;

                                    try
                                    {
                                        mass = ParseMass(parts[0], line);
                                    }
                                    catch (Exception e)
                                    {
                                        Console.WriteLine(
                                            "  Skipping spectrum with scannumber '" +
                                            currentSpectrum.ScanNumber + "'. Error in parsing mass of peak.");
                                        Console.WriteLine(e.ToString());
                                        isCorrect = false;
                                        continue;
                                    }

                                    try
                                    {
                                        intensity = ParseIntensity(parts[1], line);
                                    }
                                    catch (Exception e)
                                    {
                                        Console.WriteLine(
                                            "  Skipping spectrum with scannumber '" +
                                            currentSpectrum.ScanNumber + "'. Error in parsing intensity of peak.");
                                        Console.WriteLine(e.ToString());
                                        isCorrect = false;
                                        continue;
                                    }

                                    int key = MSAMANDA_CHEMICALUTILS.ChemicalUtils.GetMassIndex(mass);
                                    if (peaks.ContainsKey(key))
                                    {
                                        if (peaks[key] < intensity)
                                            peaks[key] = intensity;
                                    }
                                    else
                                    {
                                        peaks.Add(key, intensity);
                                    }

                                    currentSpectrum.FragmentsPeaks.Add(
                                        GenerateFragmentPeak(mass, peaks[key], charge));
                                }
                                else
                                {
                                    Console.WriteLine(
                                        "  Skipping spectrum with scannumber '" +
                                        currentSpectrum.ScanNumber + "'. Error in parsing peak.");
                                    isCorrect = false;
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error parsing mgf file at or after spectrum '" + title + "'.");
                Console.WriteLine(ex.ToString());
                spectra.Clear();
                throw;
            }

            return spectra;
        }

        private static Spectrum GenerateSpectrum(List<AMassCentroid> fragmentsPeaks, int numb, int scanId, double rt, double mOverZ, int charge)
        {
            Spectrum s = new Spectrum
            {
                FragmentsPeaks = fragmentsPeaks,
                ScanNumber = numb,
                SpectrumId = scanId,
                RT = rt,
                ImmuneMasses = new SortedSet<double>(),
                ImmunePeaks = new Dictionary<int, double>()
            };
            s.Precursor.SetMassCharge(mOverZ, charge, true);
            return s;
        }

        public static AMassCentroid GenerateFragmentPeak(double position, double intensity, int charge)
        {
            AMassCentroid amass = new AMassCentroid
            {
                Position = position,
                Intensity = intensity,
                Charge = (short)charge
            };

            return amass;
        }

        private static double ParseIntensity(string replace, string line)
        {
            string ReplaceDecimalSeperator(string toReplace)
            {
                return toReplace.Replace(",", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal)
                    .Replace(".", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal);
            }

            string i = ReplaceDecimalSeperator(replace);
            double intensity = Double.Parse(i);
            return intensity;
        }

        private static double ParseMass(string replace, string line)
        {
            string ReplaceDecimalSeperator(string toReplace)
            {
                return toReplace.Replace(",", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal)
                    .Replace(".", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal);
            }

            string m = ReplaceDecimalSeperator(replace);
            double mass = Double.Parse(m);
            return mass;
        }

        private static int ParseCharge(string line)
        {
            int inx = line.IndexOf("+", StringComparison.Ordinal);
            string m = line.Substring(7);
            if (inx != -1)
                m = line.Substring(7, inx - 7);
            if (!String.IsNullOrEmpty(m))
            {
                return Int32.Parse(m);
            }

            return 0;
        }

        private static double ParseRt(string line)
        {
            string ReplaceDecimalSeperator(string toReplace)
            {
                return toReplace.Replace(",", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal)
                    .Replace(".", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal);
            }

            int inx = line.IndexOf("=", StringComparison.Ordinal);
            string m = line.Substring(inx + 1);
            m = ReplaceDecimalSeperator(m);

            bool ok = Double.TryParse(m, out var rt);
            if (!ok) throw new Exception("Error in parsing: " + line);

            return rt;
        }

        private static int ParseScanNumber(string line)
        {
            int inx = line.IndexOf("=", StringComparison.Ordinal);
            string m = line.Substring(inx + 1);
            if (m.Contains("MSMS:", StringComparison.Ordinal))
            {
                int idxMSMS = m.IndexOf("MSMS:", StringComparison.Ordinal);
                m = m.Substring(idxMSMS + 5);
            }

            return Int32.Parse(m.Trim());
        }

        private static double ParseMOverZ(string line)
        {
            string ReplaceDecimalSeperator(string toReplace)
            {
                return toReplace.Replace(",", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal)
                    .Replace(".", NumberFormatInfo.CurrentInfo.NumberDecimalSeparator, StringComparison.Ordinal);
            }

            double GetDoubleFromString(string p)
            {
                if (string.IsNullOrEmpty(p)) return -1;

                try
                {
                    return p.Contains(".", StringComparison.Ordinal)
                        ? double.Parse(p, CultureInfo.InvariantCulture)
                        : double.Parse(ReplaceDecimalSeperator(p));
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.ToString());
                    return -1;
                }
            }

            int inx = line.IndexOf(" ", StringComparison.Ordinal);
            if (inx == -1)
            {
                inx = line.IndexOf("\t", StringComparison.Ordinal);
            }
            var m = (inx == -1) ? line.Substring(8) : line.Substring(8, inx - 8);

            m = ReplaceDecimalSeperator(m);
            var mOverZ = GetDoubleFromString(m);
            return mOverZ;
        }

        private static bool LineCanBeIgnored(string line)
        {
            // line is not ignored 
            //BEGIN IONS, END IONS, TITLE, SCANS, RTINSECONDS, PEPMASS , numbers
            if (string.IsNullOrWhiteSpace(line))
                return true;

            var text = line.TrimStart().ToUpper();
            // comment lines
            if (text.StartsWith("#", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("!", StringComparison.Ordinal))
                return true;
            if (text.StartsWith(";", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("/", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("_", StringComparison.Ordinal))
                return true;
            // unused info in lines
            if (text.StartsWith("MASS", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("INSTRUMENT", StringComparison.Ordinal))
                return true;
            // new ignored values from http://www.matrixscience.com/help/data_file_help.html#RULES
            if (text.StartsWith("ACCESSION", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("CLE", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("COM", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("CUTOUT", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("COMP", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("DB", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("DECOY", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("ERRORTOLERANT", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("ETAG", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("FORMAT", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("FRAMES", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("IT_MODS", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("ITOL", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("ITOLU", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("LIBRARY_SEARCH", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("LOCUS", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("MODS", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("MULTI_SITE_MODS", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("PEP_ISOTOPE_ERROR", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("PFA", StringComparison.Ordinal))
                return true;
            // only used in .pks, .xml
            if (text.StartsWith("PRECURSOR", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("QUANTIFICATION", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("RAWFILE", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("RAWSCANS", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("REPORT", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("REPTYPE", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("SEARCH", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("SEG", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("SEQ", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("TAG", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("TAXONOMY", StringComparison.Ordinal))
                return true;
            if (text.StartsWith("TOL", StringComparison.Ordinal))
                return true;
            // combines USER00 - USER12, USEREMAIL, USERNAME
            if (text.StartsWith("USER", StringComparison.Ordinal))
                return true;

            return false;
        }
    }
}
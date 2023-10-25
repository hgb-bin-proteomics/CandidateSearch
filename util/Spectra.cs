namespace CandidateSearch.util
{
    public class Spectrum
    {
        public double[] mz;
        public double[] intensity;
        public int scanNumber;

        public Spectrum(double[] Mz, double[] Intensity, int ScanNumber) 
        {
            mz = Mz;
            intensity = Intensity;
            scanNumber = ScanNumber;
            Array.Sort(mz, intensity);
        }

        public int[] getEncoding(int massRange = 1300, int massMultiplier = 100)
        {
            var encoding = new List<int>();

            for (int i = 0; i < mz.Length; i++)
            {
                if (mz[i] < massRange)
                {
                    encoding.Add((int) (mz[i] * massMultiplier));
                }
            }

            return encoding.ToArray();
        }
    }

    public static class MGFReader
    {
        public static List<Spectrum> readMGF(string filename, bool deconvolute = true)
        { 
            var MSAMANDA_spectra = MSAMANDA_MGFPARSER.MGFParser.ParseNextSpectra(filename);

            var spectra = new List<Spectrum>();

            foreach(var spectrum in MSAMANDA_spectra)
            {
                var mz = new List<double>();
                var intensity = new List<double>();

                foreach (var peak in spectrum.FragmentsPeaks)
                {
                    mz.Add(peak.Position);
                    intensity.Add(peak.Intensity);
                }

                spectra.Add(new Spectrum(mz.ToArray(), intensity.ToArray(), spectrum.ScanNumber));
            }

            return spectra;
        }
    }
}


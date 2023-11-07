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

        public int[] getEncoding(int massRange = 5000, int massMultiplier = 100)
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

    public class Peak
    {
        public double mz { get; set; }
        public double intensity { get; set; }
        public int charge { get; set; }
        public int envelope { get; set; }

        public Peak(double Mz, double Intensity, int Charge, int Envelope)
        {
            mz = Mz; 
            intensity = Intensity; 
            charge = Charge; 
            envelope = Envelope;
        }
    }

    public static class MGFReader
    {
        public static List<Spectrum> readMGF(string filename, Settings settings)
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

                var mzArray = mz.ToArray();
                var intensityArray = intensity.ToArray();

                if (settings.DECONVOLUTE_SPECTRA)
                {
                    SpectrumPreprocessor.deconvolute(ref mzArray, ref intensityArray, spectrum.Precursor.Charge, settings.TOLERANCE);
                }

                spectra.Add(new Spectrum(mzArray, intensityArray, spectrum.ScanNumber));
            }

            return spectra;
        }
    }
}


namespace CandidateSearch.util
{
    /// <summary>
    /// Simplified class implemention of a mass spectrum.
    /// </summary>
    public class Spectrum
    {
        /// <summary>
        /// Array containing m/z values of centroid peaks.
        /// </summary>
        public double[] mz { get; }
        /// <summary>
        /// Array containing intensities of centroid peaks.
        /// </summary>
        public double[] intensity { get; }
        /// <summary>
        /// The scan number of the spectrum.
        /// </summary>
        public int scanNumber { get; }

        /// <summary>
        /// Constructor to create a new spectrum.
        /// </summary>
        /// <param name="Mz">Array containing m/z values of centroid peaks.</param>
        /// <param name="Intensity">Array containing intensities of centroid peaks.</param>
        /// <param name="ScanNumber">The scan number of the spectrum.</param>
        public Spectrum(double[] Mz, double[] Intensity, int ScanNumber) 
        {
            mz = Mz;
            intensity = Intensity;
            scanNumber = ScanNumber;
            Array.Sort(mz, intensity);
        }

        /// <summary>
        /// Get the encoding vector of the spectrum.
        /// </summary>
        /// <param name="massRange">Maximum m/z that should be considered while encoding. Has to match the specifications of VectorSearch.</param>
        /// <param name="massMultiplier">Precision of the encoding. Has to match the specifications of VectorSearch.</param>
        /// <returns>The encoding vector as an integer array.</returns>
        public int[] getEncoding(int massRange = 5000, int massMultiplier = 100)
        {
            var encoding = new List<int>();

            for (int i = 0; i < mz.Length; i++)
            {
                if (mz[i] < massRange)
                {
                    encoding.Add((int) Math.Round(mz[i] * massMultiplier));
                }
            }

            return encoding.ToArray();
        }
    }

    /// <summary>
    /// Reader class to read mgf files.
    /// </summary>
    public static class MGFReader
    {
        /// <summary>
        /// Reads the specified mgf file and returns a list of spectra.
        /// </summary>
        /// <param name="filename">The name of the mgf file.</param>
        /// <returns>The list of spectra read from the mgf file.</returns>
        public static List<Spectrum> readMGF(string filename)
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


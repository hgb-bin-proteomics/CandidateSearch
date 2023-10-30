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

                if (settings.DECONVOLUTE_SPECTRA)
                {
                    SpectrumPreprocessor.deconvolute(ref mz, ref intensity, spectrum.Precursor.Charge, settings.TOLERANCE);
                }

                spectra.Add(new Spectrum(mz.ToArray(), intensity.ToArray(), spectrum.ScanNumber));
            }

            return spectra;
        }
    }

    public static class SpectrumPreprocessor
    {
        public static void deconvolute(ref List<Double> mzArray, 
                                       ref List<Double> intensityArray, 
                                       int precursorCharge, 
                                       double tolerance,
                                       string method = "sage")
        {
            if (method == "ms_andrea")
            {
                deconvoluteMSAndrea(ref mzArray, ref intensityArray);
                return;
            }

            deconvoluteSage(ref mzArray, ref intensityArray, precursorCharge, tolerance);
            return;
        }

        public static void deconvoluteMSAndrea(ref List<Double> mzArray, 
                                               ref List<Double> intensityArray)
        {
            // todo
            throw new NotImplementedException("MS Andrea deconvolution not yet implemented! Use Sage deconvolution instead!");
        }

        public static void deconvoluteSage(ref List<Double> mzArray, 
                                           ref List<Double> intensityArray, 
                                           int precursorCharge, 
                                           double tolerance)
        {
            double PROTON = 1.007276466812;
            double NEUTRON = 1.00335;
            var peaks = new Dictionary<int, Peak>();

            for (int i = 0; i < mzArray.Count; i++)
            {
                peaks.Add(i, new Peak(mzArray[i], intensityArray[i], 0, 0));
            }

            for (int i = mzArray.Count - 1; i >= 0; i--)
            {
                int j = i.saturatedSubtraction(1);
                while (mzArray[i] - mzArray[j] <= NEUTRON + tolerance)
                {
                    var delta = mzArray[i] - mzArray[j];
                    for (int charge = 1; charge <= precursorCharge; charge++)
                    {
                        double iso = NEUTRON / (double) charge;
                        if (Math.Abs(delta - iso) <= tolerance && intensityArray[i] < intensityArray[j])
                        {
                            if (peaks[i].charge != charge)
                            {
                                continue;
                            }

                            peaks[j].intensity += peaks[i].intensity;
                            peaks[j].charge = charge;
                            peaks[i].charge = charge;
                            peaks[i].envelope = j;
                        }
                    }
                    j = j.saturatedSubtraction(1);
                    if (j == 0)
                    {
                        break;
                    }
                }
            }

            var newMzArray = new List<double>();
            var newIntensityArray = new List<double>();

            int currentEnvelope = -1;
            for (int i = 0; i < mzArray.Count; i++)
            {
                if (peaks[i].charge != 0 && peaks[i].envelope != currentEnvelope)
                {
                    newMzArray.Add((peaks[i].mz - peaks[i].charge * PROTON) * peaks[i].charge + PROTON);
                    newIntensityArray.Add(peaks[i].intensity);
                    currentEnvelope = peaks[i].envelope;
                }

            }

            mzArray = newMzArray;
            intensityArray = newIntensityArray;
        }

        public static int saturatedSubtraction(this int x, int y)
        {
            return x - y < 0 ? 0 : x - y;
        }
    }
}


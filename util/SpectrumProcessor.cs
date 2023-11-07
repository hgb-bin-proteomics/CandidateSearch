namespace CandidateSearch.util
{
    public static class SpectrumPreprocessor
    {
        public const double PROTON = 1.007276466812;
        public const double NEUTRON = 1.00335;

        public static void deconvolute(ref double[] mzArray,
                                       ref double[] intensityArray,
                                       int precursorCharge,
                                       double tolerance,
                                       string method = "sage")
        {
            if (method == "ms_andrea")
            {
                deconvoluteMSAndrea(ref mzArray, ref intensityArray, precursorCharge);
                return;
            }

            if (method == "sage")
            {
                deconvoluteSage(ref mzArray, ref intensityArray, precursorCharge, tolerance);
                return;
            }

            deconvolute(ref mzArray, ref intensityArray, precursorCharge, tolerance);
            return;
        }

        public static void deconvolute(ref double[] mzArray,
                                       ref double[] intensityArray,
                                       int precursorCharge,
                                       double tolerance)
        {
            // todo
        }

        public static void deconvoluteMSAndrea(ref double[] mzArray,
                                               ref double[] intensityArray,
                                               int precursorCharge)
        {
            MSANDREA_DECONVOLUTION.Deconvolutor.deconvolute(ref mzArray, ref intensityArray, precursorCharge);
            return;
        }

        public static void deconvoluteSage(ref double[] mzArray,
                                           ref double[] intensityArray,
                                           int precursorCharge,
                                           double tolerance)
        {
            var peaks = new Dictionary<int, Peak>();

            Array.Sort(mzArray, intensityArray);

            for (int i = 0; i < mzArray.Length; i++)
            {
                peaks.Add(i, new Peak(mzArray[i], intensityArray[i], 0, 0));
            }

            for (int i = mzArray.Length - 1; i >= 0; i--)
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
                            if (peaks[i].charge != 0)
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
            for (int i = 0; i < mzArray.Length; i++)
            {
                if (peaks[i].charge != 0 && peaks[i].envelope != currentEnvelope)
                {
                    newMzArray.Add(calculateUnchargedMass(peaks[i].mz, peaks[i].charge) + PROTON);
                    newIntensityArray.Add(peaks[i].intensity);
                    currentEnvelope = peaks[i].envelope;
                }

            }

            mzArray = newMzArray.ToArray();
            intensityArray = newIntensityArray.ToArray();

            return;
        }

        public static int saturatedSubtraction(this int x, int y)
        {
            return x - y < 0 ? 0 : x - y;
        }

        public static double calculateUnchargedMass(double mz, int charge)
        {
            return ((mz * (double) charge) - ((double) charge * PROTON));
        }
    }
}

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
                                       string method = "ms_andrea")
        {
            if (method == "ms_andrea")
            {
                deconvoluteMSAndrea(ref mzArray, ref intensityArray, precursorCharge);
                return;
            }

            return;
        }

        public static void deconvoluteMSAndrea(ref double[] mzArray,
                                               ref double[] intensityArray,
                                               int precursorCharge)
        {
            MSANDREA_DECONVOLUTION.Deconvolutor.deconvolute(ref mzArray, ref intensityArray, precursorCharge);
            return;
        }
    }
}

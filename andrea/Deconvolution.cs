using MSAMANDA_CHEMICALUTILS;

namespace MSANDREA_DECONVOLUTION
{
    public static class Deconvolutor
    {
        public static void deconvolute(ref double[] mzArray, ref double[] intensityArray, int charge)
        {
            var masses = new SortedSet<double>();
            var peaks = new Dictionary<int, double>();
            for (int i = 0; i < mzArray.Length; i++)
            {
                var key = ChemicalUtils.GetMassIndex(mzArray[i]);
                if (peaks.ContainsKey(key))
                {
                    peaks[key] += intensityArray[i];
                }
                else
                {
                    peaks.Add(key, intensityArray[i]);
                    masses.Add(mzArray[i]);
                }
            }

            deisotope(ref masses, ref peaks, charge, out SortedSet<double> doublyChargedPeaks, out SortedSet<double> triplyChargedPeaks);
            chargeDeconvolute(ref masses, ref peaks, charge, ref doublyChargedPeaks, ref triplyChargedPeaks);

            mzArray = masses.ToArray();
            intensityArray = peaks.Values.ToArray();

            return;
        }

        public static void deisotope(ref SortedSet<double> masses, 
                                     ref Dictionary<int, double> peaks, 
                                     int charge,
                                     out SortedSet<double> doublyChargedPeaks, 
                                     out SortedSet<double> triplyChargedPeaks)
        {
            doublyChargedPeaks = new SortedSet<double>();
            triplyChargedPeaks = new SortedSet<double>();
            
            double last = 0.0;
            double lastIntensity = 0.0;
            double summedIntensity = 0.0;
            int lastMassInt = 0;
            bool foundbefore = false;
            List<double> isotopicallyReducedMasses = new List<double>();
            Tolerance tol = new Tolerance(0.005, MassUnit.DA);
            var processedMasses = new List<double>();

            // Get the starting (charge) point
            int start = charge > 3 ? 3 : charge;
            for (int j = start; j > 0; --j)
            {
                Dictionary<int, double> tempPeaks = new Dictionary<int, double>(peaks.Keys.Count);
                SortedSet<double> tempMasses = new SortedSet<double>();
                SortedSet<double> tempDoublyChargedPeaks = new SortedSet<double>();
                SortedSet<double> tempTriplyChargedPeaks = new SortedSet<double>();
                last = 0;
                lastMassInt = 0;
                summedIntensity = 0.0;

                // Search for peaks displaced by the mass of an isotope (according to the charge (j) of the current iteration of the outer most loop)
                foreach (double mass in masses)
                {
                    bool alreadyProcessed = processedMasses.Contains(mass);
                    if (!alreadyProcessed)
                    {
                        bool lastAlreadyProcessed = processedMasses.Contains(last);
                        if (last > 0.0 && !lastAlreadyProcessed)
                        {
                            double c13 = ChemicalUtils.Chemicals["i"].MonoMass / (j * 1.0);
                            double lower = tol.CalculateLowerBound(mass);
                            double upper = tol.CalculateUpperBound(mass);
                            double possibleIsotope = last + c13;

                            if (possibleIsotope > lower && possibleIsotope < upper)
                            { // If the possible isotope lies within the tolerance

                                // Find the entire isotopic distribution 
                                List<double> isotopicDistribution = GetIsotopicDistribution(last, j, ref masses, ref peaks, ref processedMasses);
                                if (isotopicDistribution.Any())
                                {
                                    summedIntensity = GetSummedIntensityIsoTopicDistribution(isotopicDistribution, ref peaks, ref processedMasses);
                                    processedMasses.Remove(last); // 'last' has to be removed from the list, otherwise it will not be in the final list after the second round if the loop looking for singly charged distributions
                                    lastMassInt = ChemicalUtils.GetMassIndex(last);
                                    tempPeaks.Add(lastMassInt, summedIntensity);
                                    tempMasses.Add(last);
                                    // Doubly  and triply charged peaks used for charge reduction 
                                    if (j == 2)
                                    {
                                        tempDoublyChargedPeaks.Add(last);
                                    }
                                    else if (j == 3)
                                    {
                                        tempTriplyChargedPeaks.Add(last);
                                    }
                                }
                                else
                                { // If no isotopic distribution is found
                                    lastMassInt = ChemicalUtils.GetMassIndex(last);
                                    lastIntensity = peaks[lastMassInt];

                                    tempPeaks.Add(lastMassInt, lastIntensity);
                                    tempMasses.Add(last);
                                }
                            }
                            else
                            { // if a possible isotope for 'last' is not found 
                                lastMassInt = ChemicalUtils.GetMassIndex(last);
                                lastIntensity = peaks[lastMassInt];

                                tempPeaks.Add(lastMassInt, lastIntensity);
                                tempMasses.Add(last);
                            }
                        }
                        last = mass;
                    }
                }

                // Setting objects after the last mass 
                int lastKey = ChemicalUtils.GetMassIndex(last);
                lastIntensity = peaks[lastKey];
                tempPeaks.Add(lastKey, lastIntensity);
                tempMasses.Add(last);
                masses = tempMasses;
                peaks = tempPeaks;
                if (j == 2)
                {
                    doublyChargedPeaks = tempDoublyChargedPeaks;
                }
                else if (j == 3)
                {
                    triplyChargedPeaks = tempTriplyChargedPeaks;
                }
            }
        }

        public static List<double> GetIsotopicDistribution(double peak, 
                                                           int charge, 
                                                           ref SortedSet<double> masses, 
                                                           ref Dictionary<int, double> peaks, 
                                                           ref List<double> processedMasses)
        {
            /*
             * Inspired from Michas method in MS Annika 
             */

            var isotopicPeaksList = new List<double>();
            var isotopicDistribution = new List<double>();
            var chargeAsDouble = (double) charge;
            var tolerance = 0.005;
            var isotopeDiff = ChemicalUtils.Chemicals["i"].MonoMass; // The mass of one neutron

            // maximum intensity difference (in m/z direction) between peaks of the same isotopic envelope. Given as 100 - d * 100 % of previous peak (From Micha)
            var d = 0.2;

            var numOfIsotopicPeaks = 0;
            var intMassCurrent = 0;
            var intensityCurrent = 0.0;
            var intMassPrev = 0;
            var intensityPrev = 0.0;
            var primaryPeakInIsotopicDistribution = false;

            // If the current peak have already been reduced to the peak right before the current one. Then the empty distribution will be returned  
            if (processedMasses.Contains(peak))
            {
                return isotopicDistribution;
            }

            // Only selecting the masses that are in close proximity of the peak
            var massesSubset = masses.Where(m => m > peak && m < peak + 8.0).ToList();

            // Setting the number of isotopic peaks to look for 
            if (massesSubset.Count >= 7)
            {
                numOfIsotopicPeaks = 7;
            }
            else
            {
                numOfIsotopicPeaks = massesSubset.Count;
            }

            // Adding the peak to the list of isotopic peaks 
            isotopicPeaksList.Add(peak);

            // Looking for isotopic peaks 
            for (int i = 1; i <= numOfIsotopicPeaks; i++)
            {
                if (isotopicPeaksList.Count == i)
                { // This (hopefully) insures that the isotopic distribution has no 'gaps' 
                    foreach (double mass in massesSubset)
                    {
                        if (mass <= (peak + i * (isotopeDiff / chargeAsDouble + tolerance)) &&
                            mass >= (peak + i * (isotopeDiff / chargeAsDouble - tolerance)))
                        {
                            isotopicPeaksList.Add(mass);
                        }
                    }
                }
                else break;
            }

            while (!primaryPeakInIsotopicDistribution)
            {
                // Only if the isotopic peaks are found, post processing is performed 
                if (isotopicPeaksList.Count > 1)
                {
                    // Post processing: Checks if the peaks follow a correct distribution 
                    int indexOfMostIntense = GetIndexOfMostIntensePeak(isotopicPeaksList, peaks);

                    // If the first peak is the most intense one 
                    if (indexOfMostIntense == 0)
                    {
                        isotopicDistribution.Add(isotopicPeaksList[0]);

                        for (int i = 1; i < isotopicPeaksList.Count; i++)
                        {
                            // Checking that the next peak in the list have a 'reasonable' intensity 
                            intMassCurrent = ChemicalUtils.GetMassIndex(isotopicPeaksList[i]);
                            intensityCurrent = peaks[intMassCurrent];

                            intMassPrev = ChemicalUtils.GetMassIndex(isotopicPeaksList[i - 1]);
                            intensityPrev = peaks[intMassPrev];
                            if (intensityCurrent < intensityPrev)
                            {
                                isotopicDistribution.Add(isotopicPeaksList[i]);
                            }
                            else break;
                        }
                    }
                    // If the first peak is not the most intense one and the detected isotopic distribution had more than two peaks and the precursor charge 
                    else if (indexOfMostIntense > 0 && isotopicPeaksList.Count > 2 && charge > 1)
                    {
                        isotopicDistribution.Add(isotopicPeaksList[indexOfMostIntense]);

                        // Checking if the peak intensities follow a correct distribution on the RIGHT side of the most intense peak
                        for (int i = indexOfMostIntense + 1; i < isotopicPeaksList.Count; i++)
                        {
                            intMassCurrent = ChemicalUtils.GetMassIndex(isotopicPeaksList[i]);
                            intensityCurrent = peaks[intMassCurrent];

                            intMassPrev = ChemicalUtils.GetMassIndex(isotopicPeaksList[i - 1]);
                            intensityPrev = peaks[intMassPrev];

                            if (intensityCurrent < intensityPrev &&
                                intensityCurrent > intensityPrev * d)
                            {
                                isotopicDistribution.Add(isotopicPeaksList[i]);
                            }
                            else break;
                        }

                        // Checking if the peak intensities follow a correct distribution on the LEFT side of the most intense peak
                        for (int i = indexOfMostIntense; i > 0; i--)
                        {
                            intMassCurrent = ChemicalUtils.GetMassIndex(isotopicPeaksList[i]);
                            intensityCurrent = peaks[intMassCurrent];

                            intMassPrev = ChemicalUtils.GetMassIndex(isotopicPeaksList[i - 1]);
                            intensityPrev = peaks[intMassPrev];

                            if (intensityCurrent > intensityPrev &&
                                intensityPrev > intensityCurrent * d)
                            {
                                isotopicDistribution.Add(isotopicPeaksList[i - 1]);
                            }
                            else
                                break;
                        }

                        // Ordering the list by increasing mass 
                        isotopicDistribution = isotopicDistribution.OrderBy(m => m).ToList();
                    }

                    // If the identified isotopic distribution does not contain the peak for which the method was called 
                    // (due to the intensities not following a correct distribution) the isotopic distribution is cleared 
                    if (!isotopicDistribution.Contains(peak))
                    {
                        isotopicDistribution.Clear();
                        isotopicPeaksList.RemoveRange(indexOfMostIntense, isotopicPeaksList.Count - indexOfMostIntense);
                    }
                }
                else
                {
                    isotopicDistribution.Add(peak);
                }

                primaryPeakInIsotopicDistribution = isotopicDistribution.Contains(peak);

            }

            return isotopicDistribution;
        }

        public static int GetIndexOfMostIntensePeak(List<double> isotopicPeaksList, 
                                                    Dictionary<int, double> peaks)
        {
            /*
             * Inspired by Michas method in MS Annika 
             */

            int maxIndex = 0;
            double maxIntensity = 0.0;

            for (int i = 0; i < isotopicPeaksList.Count; i++)
            {
                int intMass = ChemicalUtils.GetMassIndex(isotopicPeaksList[i]);
                if (peaks[intMass] > maxIntensity)
                {
                    maxIndex = i;
                    maxIntensity = peaks[intMass];
                }
            }

            return maxIndex;
        }

        public static double GetSummedIntensityIsoTopicDistribution(List<double> isotopicDistribution, 
                                                                    ref Dictionary<int, double> peaks, 
                                                                    ref List<double> processedMasses)
        {
            double intensitySum = 0.0;

            foreach (double mass in isotopicDistribution)
            {
                var intMass = ChemicalUtils.GetMassIndex(mass);
                intensitySum += peaks[intMass];
                processedMasses.Add(mass); // Moved here from the 'GetIsotopicDistribution' method
            }

            return intensitySum;
        }

        public static void chargeDeconvolute(ref SortedSet<double> masses, 
                                             ref Dictionary<int, double> peaks, 
                                             int charge,
                                             ref SortedSet<double> doubleChargedMasses, 
                                             ref SortedSet<double> triplyChargedMasses)
        {
            //SortedSet<double> masses = currentSpectrum.Masses2;
            var massesNew = new SortedSet<double>();
            var massesSubset = new SortedSet<double>();

            var peaksNew = new Dictionary<int, double>();
            var peaksSubset = new Dictionary<int, double>();
       
            var protonMassDouble = ChemicalUtils.Chemicals["p"].MonoMass;

            if (charge > 1)
            {
                // Preparing the new 'peaks' and 'masses' 
                foreach (var mass in masses)
                {
                    massesNew.Add(mass);
                    int massAsInt = ChemicalUtils.GetMassIndex(mass);
                    peaksNew.Add(massAsInt, peaks[massAsInt]);
                }

                // Dealing with doubly charged peaks 
                if (doubleChargedMasses != null)
                {
                    foreach (double doubleChargedMass in doubleChargedMasses)
                    {
                        int intMass = ChemicalUtils.GetMassIndex(doubleChargedMass);
                        if (peaksNew.ContainsKey(intMass))
                        {
                            double doubleSingleChargeMass = (doubleChargedMass * 2.0) - protonMassDouble;
                            int intSingleChargeMass = ChemicalUtils.GetMassIndex(doubleSingleChargeMass);

                            // Dealing with peaks dictionary 
                            if (peaks.ContainsKey(intSingleChargeMass))
                            {
                                double summedIntensity = peaks[intMass] + peaks[intSingleChargeMass];

                                List<double> singlyChargedPeak = new List<double>();

                                // Finding the singly charged peak
                                double tolerance = 0.005;
                                while (singlyChargedPeak.Count != 1)
                                {
                                    singlyChargedPeak = masses.Where(p => p - tolerance <= doubleSingleChargeMass && p + tolerance >= doubleSingleChargeMass).ToList();
                                    tolerance = -0.001;
                                }

                                // Adding the singly charged peak with the summed mass to the peaks dictionary 
                                if (singlyChargedPeak.Count == 1)
                                {
                                    doubleSingleChargeMass = singlyChargedPeak[0];
                                    intSingleChargeMass = ChemicalUtils.GetMassIndex(singlyChargedPeak[0]);
                                    peaksSubset.Add(intSingleChargeMass, summedIntensity);
                                }


                            }
                            else if (peaksSubset.ContainsKey(intSingleChargeMass))
                            {
                                double summedIntensity = peaks[intMass] + peaksSubset[intSingleChargeMass];
                                peaksSubset[intSingleChargeMass] = summedIntensity;
                            }
                            else
                            {
                                peaksSubset.Add(intSingleChargeMass, peaks[intMass]);
                            }
                            peaksNew.Remove(intMass);

                            // Dealing with masses list  
                            massesSubset.Add(doubleSingleChargeMass);
                            massesNew.Remove(doubleChargedMass);
                        }
                    }
                }

                //Dealing with triply charged peaks triply charged peaks
                if (triplyChargedMasses != null)
                {
                    foreach (double triplyChargedMass in triplyChargedMasses)
                    {
                        int intMass = ChemicalUtils.GetMassIndex(triplyChargedMass);
                        if (peaksNew.ContainsKey(intMass))
                        {
                            // Dealing with peaks dictionary 
                            double doubleSingleChargeMass = (triplyChargedMass * 3.0) - (protonMassDouble * 2.0);
                            int intSingleChargeMass = ChemicalUtils.GetMassIndex(doubleSingleChargeMass);

                            if (peaks.ContainsKey(intSingleChargeMass))
                            {
                                double summedIntensity = peaks[intMass] + peaks[intSingleChargeMass];

                                List<double> singlyChargedPeak = new List<double>();

                                double tolerance = 0.005;

                                while (singlyChargedPeak.Count != 1)
                                {
                                    singlyChargedPeak = masses.Where(p => p - tolerance <= doubleSingleChargeMass && p + tolerance >= doubleSingleChargeMass).ToList();
                                    tolerance = -0.001;
                                }

                                peaksSubset[intSingleChargeMass] = summedIntensity;
                            }
                            else if (peaksSubset.ContainsKey(intSingleChargeMass))
                            {
                                double summedIntensity = peaks[intMass] + peaksSubset[intSingleChargeMass];
                                peaksSubset[intSingleChargeMass] = summedIntensity;
                            }
                            else
                            {
                                peaksSubset.Add(intSingleChargeMass, peaksNew[intMass]);
                            }
                            peaksNew.Remove(intMass);

                            // Dealing with masses list
                            massesSubset.Add(doubleSingleChargeMass);
                            massesNew.Remove(triplyChargedMass);
                        }
                    }
                }


                // Adding the singly charged ions to 'peaks' and 'masses'
                if (peaksSubset.Count == massesSubset.Count)
                {
                    // peaks dictionary
                    foreach (var kvp in peaksSubset)
                    {
                        if (peaksNew.ContainsKey(kvp.Key))
                        {
                            peaksNew[kvp.Key] = kvp.Value;
                        }
                        else
                        {
                            peaksNew.Add(kvp.Key, kvp.Value);
                        }
                    }

                    // Masses list 
                    foreach (var mass in massesSubset)
                    {
                        massesNew.Add(mass);
                    }
                }
            }
   
            masses = massesNew;
            peaks = peaksNew;
        }
    }

    public enum MassUnit
    {
        DA,
        PPM,
        MMU
    };

    public class Tolerance
    {
        public double Value { get; private set; }
        public MassUnit Unit { get; private set; }
        double ppmLowerBoundFactor;
        double ppmUpperBoundFactor;
        const double mmuFactor = 1000.0;

        public Tolerance(double val, MassUnit un)
        {

            Unit = un;
            Value = val;

            const double ppmFactor = 1.0 / 1000000.0;

            if (un == MassUnit.PPM)
            {
                ppmLowerBoundFactor = 1.0 / (1.0 + Value * ppmFactor);
                ppmUpperBoundFactor = 1.0 / (1.0 - Value * ppmFactor);
            }
            else
            {
                ppmLowerBoundFactor = 1.0;
                ppmUpperBoundFactor = 1.0;
            }
        }

        public double CalculateLowerBound(double mass)
        {

            switch (Unit)
            {
                case MassUnit.DA:
                    return (mass - Value);
                case MassUnit.PPM:
                    return (mass * ppmLowerBoundFactor);
                case MassUnit.MMU:
                    return (mass - Value / mmuFactor);
            }

            return 0.0;
        }

        public double CalculateUpperBound(double mass)
        {

            switch (Unit)
            {
                case MassUnit.DA:
                    return (mass + Value);
                case MassUnit.PPM:
                    return (mass * ppmUpperBoundFactor);
                case MassUnit.MMU:
                    return (mass + Value / mmuFactor);
            }

            return 0.0;
        }
    }
}

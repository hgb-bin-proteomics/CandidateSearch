using MSAMANDA_CHEMICALUTILS;

namespace MSANDREA_DECONVOLUTION
{
    public static class Deconvolutor
    {
        public static void deconvolute(ref double[] mzArray, ref double[] intensityArray, int charge)
        {

        }

        public static void chargeDeconvolute(ref double[] mzArray, ref double[] intensityArray, int charge,
                                             ref SortedSet<double> doubleChargedMasses, ref SortedSet<double> triplyChargedMasses)
        {
            //SortedSet<double> masses = currentSpectrum.Masses2;
            var massesNew = new SortedSet<double>();
            var massesSubset = new SortedSet<double>();

            //Dictionary<int, double> peaks = currentSpectrum.Peaks;
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
   
            mzArray = massesNew.ToArray();
            intensityArray = peaksNew.Values.ToArray();
        }
    }
}

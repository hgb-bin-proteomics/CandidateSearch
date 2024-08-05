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

using System.Text;
using MSAMANDA_CHEMICALUTILS;

namespace MSAMANDA_IONCALCULATION
{
    public static class IonCalculator
    {
        static double massH_MonoMass = ChemicalUtils.Chemicals["H"].MonoMass;
        static double massProton_MonoMass = ChemicalUtils.Chemicals["p"].MonoMass;
        static double massO_MonoMass = ChemicalUtils.Chemicals["O"].MonoMass;
        static double massC_MonoMass = ChemicalUtils.Chemicals["C"].MonoMass;
        static double massN_MonoMass = ChemicalUtils.Chemicals["N"].MonoMass;

        static double massH_AvgMass = ChemicalUtils.Chemicals["H"].AvgMass;
        static double massProton_AvgMass = ChemicalUtils.Chemicals["p"].AvgMass;
        static double massO_AvgMass = ChemicalUtils.Chemicals["O"].AvgMass;
        static double massC_AvgMass = ChemicalUtils.Chemicals["C"].AvgMass;
        static double massN_AvgMass = ChemicalUtils.Chemicals["N"].AvgMass;

        public static bool CalculateIons(out double[] outIonsNoNL, 
                                         out IonWithNL[] outIonsWithNL,
                                         byte[] sequence,
                                         double mass,
                                         int charge, 
                                         Modification[] mods, 
                                         int maxNumberNeutralLoss, 
                                         int maxNumberNeutralLossModifications,
                                         double lowerBound,
                                         double upperBound, 
                                         bool mono, 
                                         string maxAllowedChargeState)
        {
            var setting = new InstrumentSetting() { UseBIon = true, UseYIon = true };

            var AminoAcids = new AminoAcid[byte.MaxValue];
            foreach (AminoAcid aa in ChemicalUtils.AminoAcids.Values)
            {
                AminoAcids[(byte) aa.OneLetterCode] = aa;
            }

            var MaxNrWaterLosses = 0;
            var MaxNrAmmoniaLosses = 0;

            outIonsNoNL = null;
            outIonsWithNL = null;

            int numIonsExpected = sequence.Length;
            {
                int numExpectedFactor = 0;
                if (setting.UseAIon) ++numExpectedFactor;
                if (setting.UseBIon) ++numExpectedFactor;
                if (setting.UseCIon) ++numExpectedFactor;
                if (setting.UseXIon) ++numExpectedFactor;
                if (setting.UseYIon) ++numExpectedFactor;
                if (setting.UseZIon) ++numExpectedFactor;
                if (setting.UseZPlusHIon) ++numExpectedFactor;
                if (setting.UseZPlusTwoHIon) ++numExpectedFactor;
                // new ions
                if (setting.UseAPlusHIon) ++numExpectedFactor;
                if (setting.UseBPlusHIon) ++numExpectedFactor;
                if (setting.UseCPlusHIon) ++numExpectedFactor;
                if (setting.UseXPlusHIon) ++numExpectedFactor;
                if (setting.UseYPlusHIon) ++numExpectedFactor;
                if (setting.UseAMinusHIon) ++numExpectedFactor;
                if (setting.UseBMinusHIon) ++numExpectedFactor;
                if (setting.UseCMinusHIon) ++numExpectedFactor;
                if (setting.UseXMinusHIon) ++numExpectedFactor;
                if (setting.UseYMinusHIon) ++numExpectedFactor;
                if (setting.UseZMinusHIon) ++numExpectedFactor;
                if (numExpectedFactor < 2)
                    numExpectedFactor = 2;

                numIonsExpected *= numExpectedFactor;
            }

            List<double> ionsNoNL = new List<double>(numIonsExpected);   // ions without neutral Losses
            List<IonWithNL> ionsWithNL = (setting.UseWaterLosses || setting.UseAmmoniaLosses) ? new List<IonWithNL>(numIonsExpected) : null;   // ions without neutral Losses

            double massH = massH_MonoMass;
            double massProton = massProton_MonoMass;
            double massO = massO_MonoMass;
            double massC = massC_MonoMass;
            double massN = massN_MonoMass;
            if (!mono)
            {
                massH = massH_AvgMass;
                massProton = massProton_AvgMass;
                massO = massO_AvgMass;
                massC = massC_AvgMass;
                massN = massN_AvgMass;
            }

            try
            {
                //fwd: n-terminal ions
                //rev: c-terminal ions
                int possibleAmoniaLosses = 0;
                int possibleWaterLosses = 0;
                int possibleAmoniaLossesRev = 0;
                int possibleWaterLossesRev = 0;
                //N-Term
                double lastMassFwd = massH + massProton;
                double lastMassInternalFwd = mass - massH - massO - massH + massProton;
                double lastMassInternalRev = mass - massH - massO - massH + massProton;

                //N-Terminal modification
                if (mods[0] != null)
                {
                    lastMassFwd += mods[0].Mass(mono);
                }
                //C-Term
                double lastMassRev = massH + massO + massProton;
                //C-Terminal modification
                if (mods[mods.Length - 1] != null)
                {
                    lastMassRev += mods[mods.Length - 1].Mass(mono);
                }
                //possible neutral losses of specific modifications 
                List<double> fwdLosses = new List<double>();
                List<double> revLosses = new List<double>();


                #region ion calculation
                int j = sequence.Length - 1;
                for (int i = 0; i < sequence.Length - 1; ++i, --j)
                {
                    byte aaFwd = sequence[i];
                    if (null == AminoAcids[aaFwd])
                        return false;
                    byte aaRev = sequence[j];
                    if (null == AminoAcids[aaRev])
                        return false;
                    AminoAcid currentAAFwd = AminoAcids[aaFwd];
                    lastMassFwd += currentAAFwd.Mass(mono);

                    //AA is modified
                    if (mods[i + 1] != null)
                        lastMassFwd += mods[i + 1].Mass(mono);

                    //Modification can have neutral loss
                    if (mods[i + 1] != null && mods[i + 1].NeutralLosses.Length > 0)
                    {
                        //HashSet<double> currentFwd = new HashSet<double>();
                        //add new neutral losses to existing losses
                        //neutral losses of modifications are either fixed or not; 
                        //if not fixed, the unimod file contains two entries for neutral losses:
                        //1. the neutral loss mass; 2. the mass 0.0, to declare neutral loss is not fixed
                        //foreach (double lossMass in fwdLosses)
                        for (int m = 0; m < mods[i + 1].NeutralLosses.Length; ++m)
                        {
                            double newLoss = mods[i + 1].NeutralLosses[m];
                            int count = 0;
                            while (newLoss > 0 && fwdLosses.Contains(newLoss) && count < maxNumberNeutralLossModifications)
                            {
                                newLoss += newLoss;
                                ++count;
                            }
                            if (count < maxNumberNeutralLossModifications)
                            {
                                if (!fwdLosses.Contains(newLoss))
                                    fwdLosses.Add(newLoss);
                            }
                        }
                        fwdLosses.Sort();
                    }
                    //same for reverse ions
                    AminoAcid currentAARev = AminoAcids[aaRev];
                    lastMassRev += currentAARev.Mass(mono);

                    if (mods[j + 1] != null)
                    {
                        lastMassRev += mods[j + 1].Mass(mono);
                    }
                    if (mods[j + 1] != null && mods[j + 1].NeutralLosses.Length > 0)
                    {
                        for (int m = 0; m < mods[j + 1].NeutralLosses.Length; ++m)
                        {
                            double newLoss = mods[j + 1].NeutralLosses[m];
                            int count = 0;
                            while (newLoss > 0 && revLosses.Contains(newLoss) && count < maxNumberNeutralLossModifications)
                            {
                                newLoss += newLoss;
                                ++count;
                            }
                            if (count < maxNumberNeutralLossModifications)
                            {
                                if (!revLosses.Contains(newLoss))
                                    revLosses.Add(newLoss);
                            }
                        }
                        revLosses.Sort();
                    }

                    //if ion contains AA S,T,E or D water loss is possible
                    if ((aaFwd == 'S' || aaFwd == 'T' || aaFwd == 'E' || aaFwd == 'D'))
                    {
                        if (possibleWaterLosses < maxNumberNeutralLoss)
                        {
                            ++possibleWaterLosses;
                        }
                        ++MaxNrWaterLosses;
                    }
                    //if ion contains AA N,Q,R or K ammonia loss is possible
                    if ((aaFwd == 'N' || aaFwd == 'Q' || aaFwd == 'R' || aaFwd == 'K'))
                    {
                        if (possibleAmoniaLosses < maxNumberNeutralLoss)
                        {
                            ++possibleAmoniaLosses;
                        }
                        ++MaxNrAmmoniaLosses;
                    }
                    if (setting.UseImmoniumIons)
                    {
                        //calculate immonium ion
                        double massImmo = (double)(currentAAFwd.ImmoniumMass);
                        //if AA is modified, add modification
                        if (mods[i + 1] != null)
                        {
                            massImmo += mods[i + 1].Mass(mono);
                        }
                        if (InRange(massImmo, lowerBound, upperBound))
                        {
                            ionsNoNL.Add(massImmo);
                        }
                        //also immonium ion for last AA
                        if (j == sequence.Length - 1)
                        {
                            double massImmoRev = currentAARev.ImmoniumMass;
                            if (mods[j + 1] != null)
                            {
                                massImmoRev += mods[j + 1].Mass(mono);
                            }
                            if (InRange(massImmoRev, lowerBound, upperBound))
                            {
                                ionsNoNL.Add(massImmoRev);
                            }
                        }
                    }
                    if ((aaRev == 'S' || aaRev == 'T' || aaRev == 'E' || aaRev == 'D') && possibleWaterLossesRev < maxNumberNeutralLoss)
                    {
                        ++possibleWaterLossesRev;
                    }
                    if ((aaRev == 'N' || aaRev == 'Q' || aaRev == 'R' || aaRev == 'K') && possibleAmoniaLossesRev < maxNumberNeutralLoss)
                        ++possibleAmoniaLossesRev;

                    if (setting.UseInternalFragments)
                    {
                        if (i == 0)
                        {
                            lastMassInternalFwd -= lastMassFwd - lastMassRev;
                            if (mods[i + 1] != null)
                            {
                                lastMassInternalFwd -= mods[i + 1].Mass(mono);
                            }
                            if (mods[j + 1] != null)
                            {
                                lastMassInternalFwd -= mods[j + 1].Mass(mono);
                            }
                            if (InRange(lastMassInternalFwd, lowerBound, upperBound))
                            {
                                ionsNoNL.Add(lastMassInternalFwd);
                            }
                            CalculateAllInternalCombinations(lastMassInternalFwd, i, ionsNoNL, sequence, mods, lowerBound, upperBound, mono, AminoAcids);
                        }
                        else if (i > 0 && i < sequence.Length - 1)
                        {
                            lastMassInternalFwd -= currentAAFwd.Mass(mono);
                            if (mods[i + 1] != null)
                            {
                                lastMassInternalFwd -= mods[i + 1].Mass(mono);
                            }
                            if (InRange(lastMassInternalFwd, lowerBound, upperBound))
                            {
                                ionsNoNL.Add(lastMassInternalFwd);
                            }
                            CalculateAllInternalCombinations(lastMassInternalFwd, i, ionsNoNL, sequence, mods, lowerBound, upperBound, mono, AminoAcids);
                        }
                    }

                    if (setting.UseAIon)
                    {
                        //a ion         
                        double aIon = lastMassFwd - massC - massH - massO;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO, massH, massC,
                         massN, massProton, fwdLosses, setting, aIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                    if (setting.UseBIon)
                    {
                        //b ion
                        if (i > 0)
                        {
                            double bIon = lastMassFwd - massH;
                            CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO, massH, massC,
                             massN, massProton, fwdLosses, setting, bIon, lowerBound, upperBound, maxAllowedChargeState);
                        }
                    }

                    if (setting.UseCIon)
                    {
                        //c
                        double cIon = lastMassFwd + massN + massH * 2;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO, massH, massC,
                         massN, massProton, fwdLosses, setting, cIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                    if (setting.UseXIon)
                    {
                        //x
                        double xIon = lastMassRev + massC + massO - massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev, massO, massH, massC,
                          massN, massProton, revLosses, setting, xIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                    if (setting.UseYIon)
                    {
                        //y
                        double yIon = lastMassRev + massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev, massO, massH, massC,
                          massN, massProton, revLosses, setting, yIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                    if (setting.UseZIon)
                    {
                        //z
                        double zIon = lastMassRev - massN - massH * 2;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev, massO, massH, massC,
                          massN, massProton, revLosses, setting, zIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseZPlusHIon)
                    {
                        //z + 1
                        double zPlusHIon = lastMassRev - massN - massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev, massO, massH, massC,
                          massN, massProton, revLosses, setting, zPlusHIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseZPlusTwoHIon)
                    {
                        //z + 2
                        double zPlusTwoHIon = lastMassRev - massN;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev, massO, massH, massC,
                          massN, massProton, revLosses, setting, zPlusTwoHIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    // new ions
                    if (setting.UseAPlusHIon)
                    {
                        //a + 1
                        double aIon = lastMassFwd - massC - massO;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO,
                            massH, massC,
                            massN, massProton, fwdLosses, setting, aIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseAMinusHIon)
                    {
                        //a - 1
                        double aIon = lastMassFwd - massC - massO - 2 * massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO,
                            massH, massC,
                            massN, massProton, fwdLosses, setting, aIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseBPlusHIon)
                    {
                        //b + 1
                        if (i > 0)
                        {
                            double bIon = lastMassFwd;
                            CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO,
                                massH, massC,
                                massN, massProton, fwdLosses, setting, bIon, lowerBound, upperBound, maxAllowedChargeState);
                        }
                    }
                    if (setting.UseBMinusHIon)
                    {
                        //b - 1
                        if (i > 0)
                        {
                            double bIon = lastMassFwd - 2 * massH;
                            CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO,
                                massH, massC,
                                massN, massProton, fwdLosses, setting, bIon, lowerBound, upperBound, maxAllowedChargeState);
                        }
                    }
                    if (setting.UseCPlusHIon)
                    {
                        //c +1
                        double cIon = lastMassFwd + massN + 3 * massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO,
                            massH, massC,
                            massN, massProton, fwdLosses, setting, cIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseCMinusHIon)
                    {
                        //c -1
                        double cIon = lastMassFwd + massN + massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO,
                            massH, massC,
                            massN, massProton, fwdLosses, setting, cIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseXPlusHIon)
                    {
                        //x +1
                        double xIon = lastMassRev + massC + massO;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev,
                            massO, massH, massC,
                            massN, massProton, revLosses, setting, xIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseXMinusHIon)
                    {
                        //x -1
                        double xIon = lastMassRev + massC + massO - 2 * massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev,
                            massO, massH, massC,
                            massN, massProton, revLosses, setting, xIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                    if (setting.UseYPlusHIon)
                    {
                        //y +1 
                        double yIon = lastMassRev + 2 * massH;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev,
                            massO, massH, massC,
                            massN, massProton, revLosses, setting, yIon, lowerBound, upperBound, maxAllowedChargeState);
                    }
                    if (setting.UseYMinusHIon)
                    {
                        //y - 1 
                        double yIon = lastMassRev;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev,
                            massO, massH, massC,
                            massN, massProton, revLosses, setting, yIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                    if (setting.UseZMinusHIon)
                    {
                        //z - 1
                        double zIon = lastMassRev - massN - massH * 3;
                        CalculateIon(charge, ionsNoNL, ionsWithNL, possibleAmoniaLossesRev, possibleWaterLossesRev,
                            massO, massH, massC,
                            massN, massProton, revLosses, setting, zIon, lowerBound, upperBound, maxAllowedChargeState);
                    }

                }

                #endregion
                char lastAA = (char)sequence[sequence.Length - 1];
                if ((lastAA == 'S' || lastAA == 'T' || lastAA == 'E' || lastAA == 'D'))
                {
                    ++MaxNrWaterLosses;
                }
                if ((lastAA == 'N' || lastAA == 'Q' || lastAA == 'R' || lastAA == 'K'))
                    ++MaxNrAmmoniaLosses;

                return GetUniqueIonsSorted(out outIonsNoNL, out outIonsWithNL, ionsNoNL, ionsWithNL);

            }
            catch (Exception e)
            {
                throw new Exception("Error in calculating ions for peptide " + Encoding.ASCII.GetString(sequence));
            }

            return false;
        }

        static bool InRange(double v, double min, double max)
        {
            return ((v >= min) && (v <= max));
        }

        static bool GetUniqueIonsSorted(out double[] outIonsNoNL, out IonWithNL[] outIonsWithNL, List<double> ionsNoNL, List<IonWithNL> ionsWithNL)
        {
            outIonsNoNL = null;
            outIonsWithNL = null;

            if (ionsNoNL.Count > 0)
            {
                double[] aIonsNoNL = new double[ionsNoNL.Count];
                for (int i = 0; i < aIonsNoNL.Length; ++i)
                {
                    aIonsNoNL[i] = IonMassStabilizer.Stabilize(ionsNoNL[i]);
                }
                Array.Sort(aIonsNoNL);
                //Removing duplicates from array
                int k1 = 0, k2 = 1;
                for (; k2 < aIonsNoNL.Length; ++k2)
                {
                    if (aIonsNoNL[k1] < aIonsNoNL[k2])
                    {
                        ++k1;
                    }
                    else
                    {
                        ++k2;
                        break;
                    }
                }

                for (; k2 < aIonsNoNL.Length; ++k2)
                {
                    if (aIonsNoNL[k1] < aIonsNoNL[k2])
                    {
                        ++k1;
                        aIonsNoNL[k1] = aIonsNoNL[k2];
                    }
                }

                Array.Resize(ref aIonsNoNL, ++k1);
                outIonsNoNL = aIonsNoNL;
            }

            if (null != ionsWithNL && ionsWithNL.Count > 0)
            {
                IonWithNL[] aIonsWithNL = ionsWithNL.ToArray();
                FastSort(aIonsWithNL);
                {
                    //Merging duplicates from array
                    int k1 = 0, k2 = 1;
                    for (; k2 < aIonsWithNL.Length; ++k2)
                    {
                        if (aIonsWithNL[k1].Mz < aIonsWithNL[k2].Mz)
                        {
                            ++k1;
                        }
                        else
                        {
                            aIonsWithNL[k1].AddRange(aIonsWithNL[k2]);
                            ++k2;
                            break;
                        }
                    }

                    for (; k2 < aIonsWithNL.Length; ++k2)
                    {
                        if (aIonsWithNL[k1].Mz < aIonsWithNL[k2].Mz)
                        {
                            ++k1;
                            aIonsWithNL[k1] = aIonsWithNL[k2];
                        }
                        else
                        {
                            aIonsWithNL[k1].AddRange(aIonsWithNL[k2]);
                        }
                    }

                    Array.Resize(ref aIonsWithNL, ++k1);
                    outIonsWithNL = aIonsWithNL;
                }
            }

            if (null != outIonsNoNL || null != outIonsWithNL)
            {
                if (null == outIonsNoNL)
                {
                    outIonsNoNL = Array.Empty<double>();
                }
                return true;
            }
            return false;
        }

        private static void doQuicksort(IonWithNL[] elements, int left, int right)
        {
            int i = left, j = right;
            //middle pivot; Should avoid n^2 worst case for nearly sorted array
            IonWithNL pivot = elements[left + (right - left) / 2];

            while (i <= j)
            {
                while (elements[i].Mz < pivot.Mz)
                {
                    i++;
                }

                while (elements[j].Mz > pivot.Mz)
                {
                    j--;
                }

                if (i <= j)
                {
                    // Swap
                    //IonWithNL tmp = elements[i];
                    //elements[i] = elements[j];
                    //elements[j] = tmp;

                    (elements[i], elements[j]) = (elements[j], elements[i]);

                    i++;
                    j--;
                }
            }

            // Recursive calls
            if (left < j)
            {
                doQuicksort(elements, left, j);
            }

            if (i < right)
            {
                doQuicksort(elements, i, right);
            }
        }

        public static void FastSort(IonWithNL[] ions)
        {
            if (ions.Length > 1)
            {
                doQuicksort(ions, 0, ions.Length - 1);
            }
        }

        public static void CalculateAllInternalCombinations(double lastMassInternal, 
                                                            int i,
                                                            List<double> ionsNoNL,
                                                            byte[] seq, 
                                                            Modification[] mods,
                                                            double lowerBound, 
                                                            double upperBound, 
                                                            bool mono,
                                                            AminoAcid[] aminoAcids)
        {
            double currentMass = lastMassInternal;
            for (int k = seq.Length - 2; k > i + 2; --k)
            {
                currentMass -= aminoAcids[(char) seq[k]].Mass(mono);

                if (mods[k + 1] != null)
                {
                    currentMass -= mods[k + 1].Mass(mono);
                }
                if (InRange(currentMass, lowerBound, upperBound))
                {
                    ionsNoNL.Add(currentMass);
                }
            }
        }

        public static void CalculateIon(int charge, List<double> ionsNoNL, List<IonWithNL> ionsWithNL, int possibleAmoniaLosses,
                                        int possibleWaterLosses, double massO, double massH, double massC, double massN, double massProton,
                                        List<double> losses, InstrumentSetting setting, double ion,
                                        double lowerBound, double upperBound, string maxAllowedChargeState)
        {
            //add possible neutral losses of modifications
            for (int i = 0; i < losses.Count; ++i)
            {
                double newMass = ion - losses[i];
                if (InRange(newMass, lowerBound, upperBound))
                {
                    //add water or ammonia losses
                    if ((setting.UseWaterLosses && possibleWaterLosses > 0) || (setting.UseAmmoniaLosses && possibleAmoniaLosses > 0))
                    {
                        IonWithNL ionWithNL = new IonWithNL(newMass);
                        ionsWithNL.Add(ionWithNL);
                        CalculatePossibleNeutralLosses(newMass, possibleWaterLosses, possibleAmoniaLosses,
                          1, setting.UseWaterLosses, setting.UseAmmoniaLosses, massO, massH, massC,
                          massN, ionWithNL, lowerBound, upperBound);
                    }
                    else
                    {
                        ionsNoNL.Add(newMass);
                    }
                }
                if (charge > 2)
                {
                    CalculateHigherChargedIons(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO, massH, massC, massN, massProton, setting, newMass, lowerBound, upperBound, maxAllowedChargeState);
                }
            }

            if (losses.Count == 0)
            {
                if (InRange(ion, lowerBound, upperBound))
                {
                    //add water or ammonia losses
                    if ((setting.UseWaterLosses && possibleWaterLosses > 0) || (setting.UseAmmoniaLosses && possibleAmoniaLosses > 0))
                    {
                        IonWithNL ionWithNL = new IonWithNL(ion);
                        ionsWithNL.Add(ionWithNL);
                        CalculatePossibleNeutralLosses(ion, possibleWaterLosses, possibleAmoniaLosses,
                          1, setting.UseWaterLosses, setting.UseAmmoniaLosses, massO, massH, massC,
                          massN, ionWithNL, lowerBound, upperBound);
                    }
                    else
                    {
                        ionsNoNL.Add(ion);
                    }
                }
                if (charge > 2)
                {
                    CalculateHigherChargedIons(charge, ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses, massO, massH, massC, massN, massProton, setting, ion, lowerBound, upperBound, maxAllowedChargeState);
                }
            }
        }

        public static void CalculateHigherChargedIons(int charge, List<double> ionsNoNL, List<IonWithNL> ionsWithNL, int possibleAmoniaLosses,
                                                      int possibleWaterLosses, double massO, double massH, double massC, double massN, double massProton,
                                                      InstrumentSetting setting, double newMass, double lowerBound, double upperBound, string maxAllowedChargeState)
        {
            string allowedChargeStates = "+2";
            switch (maxAllowedChargeState)
            {
                case "+2":
                    AddHigherChargedIons(ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses,
                      massO, massH, massC, massN, massProton, setting.UseWaterLosses, setting.UseAmmoniaLosses,
                      newMass, 1, lowerBound, upperBound);
                    break;
                case "+3":
                    for (int i = 1; i < charge - 1 && i < 3; ++i)
                    {
                        AddHigherChargedIons(ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses,
                          massO, massH, massC, massN, massProton, setting.UseWaterLosses, setting.UseAmmoniaLosses,
                          newMass, i, lowerBound, upperBound);
                    }
                    break;
                case "+4":
                    for (int i = 1; i < charge - 1 && i < 4; ++i)
                    {
                        AddHigherChargedIons(ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses,
                          massO, massH, massC, massN, massProton, setting.UseWaterLosses, setting.UseAmmoniaLosses,
                          newMass, i, lowerBound, upperBound);
                    }
                    break;
                case "Precursor - 1":
                    for (int i = 1; i < charge - 1; ++i)
                    {
                        AddHigherChargedIons(ionsNoNL, ionsWithNL, possibleAmoniaLosses, possibleWaterLosses,
                          massO, massH, massC, massN, massProton, setting.UseWaterLosses, setting.UseAmmoniaLosses,
                          newMass, i, lowerBound, upperBound);
                    }
                    break;
            }
        }

        public static void AddHigherChargedIons(List<double> ionsNoNL, List<IonWithNL> ionsWithNL, int possibleAmoniaLosses,
                                                int possibleWaterLosses, double massO, double massH, double massC, double massN,
                                                double massProton, bool water, bool ammonia, double ion, int currentCharge,
                                                double lowerBound, double upperBound)
        {
            double chargedIon = ion + currentCharge * massProton;
            ++currentCharge;

            double myMass = chargedIon / currentCharge;
            if (InRange(myMass, lowerBound, upperBound))
            {
                if ((water && possibleWaterLosses > 0) || (ammonia && possibleAmoniaLosses > 0))
                {
                    IonWithNL ionWithNL = new IonWithNL(myMass);
                    ionsWithNL.Add(ionWithNL);
                    CalculatePossibleNeutralLosses(chargedIon, possibleWaterLosses, possibleAmoniaLosses,
                      currentCharge, water, ammonia, massO, massH,
                      massC, massN, ionWithNL, lowerBound, upperBound);
                }
                else
                {
                    ionsNoNL.Add(myMass);
                }
            }
        }

        public static void CalculatePossibleNeutralLosses(double ionMass, int possibleWaterLosses,
                                                          int possibleAmoniaLosses, int charge, bool water, bool ammonia, double massO, double massH,
                                                          double massC, double massN, IonWithNL ionWithNL,
                                                          double lowerBound, double upperBound)
        {
            //handle neutral losses

            double massWaterLoss = massH * 2 + massO;
            double massAmoniaLoss = massH * 3 + massN;

            int currentWaterLoss = 0;
            while (water && currentWaterLoss < possibleWaterLosses)
            {
                ++currentWaterLoss;
                double chargedIon = ionMass - currentWaterLoss * massWaterLoss;
                double myMass = IonMassStabilizer.Stabilize(chargedIon / charge);
                //if(InRange( myMass, lowerBound, upperBound ))
                ionWithNL.Add(myMass);
            }

            int currentAmoniaLoss = 0;
            while (ammonia && currentAmoniaLoss < possibleAmoniaLosses)
            {
                ++currentAmoniaLoss;
                double chargedIon = ionMass - currentAmoniaLoss * massAmoniaLoss;
                double myMass = IonMassStabilizer.Stabilize(chargedIon / charge);
                //if(InRange( myMass, lowerBound, upperBound ))
                ionWithNL.Add(myMass);
            }
        }
    }

    public class IonWithNL : List<double>
    {
        public double Mz { get; set; }


        public IonWithNL(double baseMZ)
          : base(2)
        {
            Mz = IonMassStabilizer.Stabilize(baseMZ);
        }
    }

    internal static class IonMassStabilizer
    {
        const int SignificantBitsStabilizer = 47; // Approximately 14 decimal significant digits

        const ulong RoundingAdditionStabilizer = unchecked(1uL << (52 - SignificantBitsStabilizer - 1));
        const ulong ClearBitsStabilizerMask = unchecked(0xFFFFFFFFFFFFFFFFuL << (52 - SignificantBitsStabilizer));

        /*
        Implements super-fast base2 nearest-value rounding of double precision values
        to specified number of significant 2-base digits 
        Utilizes specificity of IEEE 754 double-precision format
        Works correctly in all cases except double.Max, NaN, +/-Inf etc in input
        Specially tested correctness when 52bit mantissa is all 1 (0xFFFFFFFFFFFFF)
        By benchmarks works 7-15 times faster than Math.Round()
        */
        static internal double Stabilize(double value)
        {
            #if !DISABLE_FP_STABILIZER
            unchecked
            {
                ulong im = (RoundingAdditionStabilizer + (ulong)BitConverter.DoubleToInt64Bits(value)) & ClearBitsStabilizerMask;
                return BitConverter.Int64BitsToDouble((long)im);
            }
            #else
                return value;
            #endif
        }
    }

    public class Modification
    {
        public readonly string Title; //{ get; private set; }
        public readonly string FullName; // { get; private set; }
        private readonly double MonoMass; // { get; private set; }
        private readonly double AvgMass; // { get; private set; }
        public readonly char AA; // { get; private set; }
        public double[] NeutralLosses { get; set; }
        public bool NTerminal { get; set; }
        public bool CTerminal { get; set; }
        public bool ProteinTerminus { get; set; }
        public int MaxOccurrence { get; set; }

        public int ID { get; set; }

        public double Mass(bool mono)
        {
            return mono ? MonoMass : AvgMass;
        }

        public bool Fixed { get; set; }
        public bool SemiFixed { get; set; }

        public Modification(string title, string name, double mono, double avg, char aa, bool fix, double[] neutralLosses,
            bool nTerminal, bool cTerminal, int id, bool protein, int maxOccurrence, bool semiFixed = false)
        {
            Title = title;
            FullName = name;
            MonoMass = mono;
            AvgMass = avg;
            AA = aa;
            Fixed = fix;
            SemiFixed = semiFixed;
            NeutralLosses = neutralLosses;
            NTerminal = nTerminal;
            CTerminal = cTerminal;
            ID = id;
            ProteinTerminus = protein;
            MaxOccurrence = maxOccurrence;
        }

        public Modification(Modification m)
        {
            Title = m.Title;
            FullName = m.FullName;
            MonoMass = m.MonoMass;
            AvgMass = m.AvgMass;
            AA = m.AA;
            Fixed = m.Fixed;
            SemiFixed = m.SemiFixed;
            NeutralLosses = m.NeutralLosses;
            NTerminal = m.NTerminal;
            CTerminal = m.CTerminal;
            ID = m.ID;
            ProteinTerminus = m.ProteinTerminus;
            MaxOccurrence = m.MaxOccurrence;
        }

        public Modification(Modification m, bool fix, bool nTerminal, bool cTerminal, bool protein, int maxOccurrence)
        {
            Title = m.Title;
            FullName = m.FullName;
            MonoMass = m.MonoMass;
            AvgMass = m.AvgMass;
            AA = m.AA;
            Fixed = fix;
            SemiFixed = m.SemiFixed;
            NeutralLosses = m.NeutralLosses;
            NTerminal = nTerminal;
            CTerminal = cTerminal;
            ID = m.ID;
            ProteinTerminus = protein;
            MaxOccurrence = maxOccurrence;
        }

        public override string ToString()
        {
            string name = Title + "(";
            if (CTerminal)
            {
                name += "C-Term)";
            }
            else if (NTerminal)
            {
                name += "N-Term)";
            }
            else
            {
                name += AA + ")";
            }
            return name;
        }

        #region oldstuff
        public static string GetSaveString(Modification modif)
        {
            StringBuilder builder = new StringBuilder();
            builder.Append(modif.Title).Append("°");
            builder.Append(modif.FullName).Append("°");
            builder.Append(modif.MonoMass).Append("°");
            builder.Append(modif.AvgMass).Append("°");
            builder.Append(modif.AA).Append("°");
            builder.Append(modif.Fixed).Append("°");
            foreach (double nl in modif.NeutralLosses)
            {
                builder.Append(nl).Append(":");
            }
            builder.Append("°");
            builder.Append(modif.NTerminal).Append("°");
            builder.Append(modif.CTerminal).Append("°");
            builder.Append(modif.ID).Append("°");
            builder.Append(modif.ProteinTerminus);
            return builder.ToString();
        }

        //public static Modification GetModifOfSaveString(string savedModif) {
        //  string[] parts = savedModif.Split('°');
        //  double monoMass = double.Parse(parts[2]);
        //  double avgMass = double.Parse(parts[3]);
        //  bool fixedModif = bool.Parse(parts[5]);
        //  string[] nls = parts[6].Split(new char[] {':'}, StringSplitOptions.RemoveEmptyEntries);
        //  double[] neutralLosses = new double[nls.Length];
        //  for (int i = 0; i < nls.Length; ++i) {
        //    double elem = double.Parse(nls[i]);
        //    neutralLosses[i] = elem;
        //  }
        //  bool nterm = bool.Parse(parts[7]);
        //  bool cterm = bool.Parse(parts[8]);
        //  int id = Int32.Parse(parts[9]);
        //  bool protein = bool.Parse(parts[10]);
        //  Modification m = new Modification(parts[0], parts[1], monoMass, avgMass, char.Parse(parts[4]), fixedModif, neutralLosses, nterm, cterm, id, protein);
        //  return m;
        //}
        #endregion
    }

    public class InstrumentSetting
    {
        public InstrumentSetting()
        {
            Name = "";
            UseAIon = false;
            UseAmmoniaLosses = false;
            UseBIon = false;
            UseCIon = false;
            UseImmoniumIons = false;
            UseWaterLosses = false;
            UseXIon = false;
            UseYIon = false;
            UseZIon = false;
            UseZPlusHIon = false;
            UseZPlusTwoHIon = false;
            UseInternalFragments = false;
            // new Ions
            UseAMinusHIon = false;
            UseBPlusHIon = false;
            UseBMinusHIon = false;
            UseCPlusHIon = false;
            UseCMinusHIon = false;
            UseXPlusHIon = false;
            UseXMinusHIon = false;
            UseYPlusHIon = false;
            UseYMinusHIon = false;
            UseZMinusHIon = false;
        }

        public InstrumentSetting Clone()
        {
            InstrumentSetting instr = new InstrumentSetting
            {
                Name = this.Name,
                UseAIon = this.UseAIon,
                UseAmmoniaLosses = this.UseAmmoniaLosses,
                UseBIon = this.UseBIon,
                UseCIon = this.UseCIon,
                UseImmoniumIons = this.UseImmoniumIons,
                UseWaterLosses = this.UseWaterLosses,
                UseXIon = this.UseXIon,
                UseYIon = this.UseYIon,
                UseZIon = this.UseZIon,
                UseZPlusHIon = this.UseZPlusHIon,
                UseZPlusTwoHIon = this.UseZPlusTwoHIon,
                UseInternalFragments = this.UseInternalFragments,
                // new ions
                UseAPlusHIon = UseAPlusHIon,
                UseAMinusHIon = UseAMinusHIon,
                UseBPlusHIon = UseBPlusHIon,
                UseBMinusHIon = UseBMinusHIon,
                UseCPlusHIon = UseCPlusHIon,
                UseCMinusHIon = UseCMinusHIon,
                UseXPlusHIon = UseXPlusHIon,
                UseXMinusHIon = UseXMinusHIon,
                UseYPlusHIon = UseYPlusHIon,
                UseYMinusHIon = UseYMinusHIon,
                UseZMinusHIon = UseZMinusHIon
            };
            return instr;
        }

        public string Name { get; set; }
        public bool UseAIon { get; set; }
        public bool UseBIon { get; set; }
        public bool UseYIon { get; set; }
        public bool UseWaterLosses { get; set; }
        public bool UseAmmoniaLosses { get; set; }
        public bool UseImmoniumIons { get; set; }
        public bool UseZIon { get; set; }
        public bool UseZPlusHIon { get; set; }
        public bool UseZPlusTwoHIon { get; set; }
        public bool UseXIon { get; set; }
        public bool UseCIon { get; set; }
        public bool UseInternalFragments { get; set; }

        // new Ions
        public bool UseAPlusHIon { get; set; }
        public bool UseAMinusHIon { get; set; }
        public bool UseBPlusHIon { get; set; }
        public bool UseBMinusHIon { get; set; }
        public bool UseCPlusHIon { get; set; }
        public bool UseCMinusHIon { get; set; }
        public bool UseXPlusHIon { get; set; }
        public bool UseXMinusHIon { get; set; }
        public bool UseYPlusHIon { get; set; }
        public bool UseYMinusHIon { get; set; }
        public bool UseZMinusHIon { get; set; }

    }
}

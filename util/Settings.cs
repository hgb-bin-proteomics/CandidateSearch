﻿using System.Text;

namespace CandidateSearch.util
{
    public class Settings
    {
        // config digestion
        public int MAX_CLEAVAGES { get; set; } // default: 2
        public int MIN_PEP_LENGTH { get; set; } // default: 5
        public int MAX_PEP_LENGTH { get; set; } // default: 30

        // config ion calculation
        public int MAX_CHARGE { get; set; } // default: 3
        public int MAX_NEUTRAL_LOSSES { get; set; } // default: 1
        public int MAX_NEUTRAL_LOSS_MODS { get; set; } // default: 2
        public string MAX_ALLOWED_CHARGE { get; set; } // default: "+3";

        // config vector search
        public int TOP_N { get; set; } // default: 100
        public float TOLERANCE { get; set; } // default: 0.02f
        public bool NORMALIZE { get; set; } // default: true
        public bool USE_GAUSSIAN { get; set; } // default: true

        public Settings(int MaxCleavages = 2, int MinPepLength = 5, int MaxPepLength = 30, 
                        int MaxCharge = 3, int MaxNeutralLosses = 1, int MaxNeutralLossMods = 2, string MaxAllowedCharge = "+3",
                        int TopN = 100, float Tolerance = 0.02f, bool Normalize = true, bool UseGaussian = true)
        {
            MAX_CLEAVAGES = MaxCleavages;
            MIN_PEP_LENGTH = MinPepLength;
            MAX_PEP_LENGTH = MaxPepLength;

            MAX_CHARGE = MaxCharge;
            MAX_NEUTRAL_LOSSES = MaxNeutralLosses;
            MAX_NEUTRAL_LOSS_MODS = MaxNeutralLossMods;
            MAX_ALLOWED_CHARGE = MaxAllowedCharge;

            TOP_N = TopN;
            TOLERANCE = Tolerance;
            NORMALIZE = Normalize;
            USE_GAUSSIAN = UseGaussian;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("---- SETTINGS: ----\n");
            sb.Append($"MAX_CLEAVAGES: {MAX_CLEAVAGES}\n");
            sb.Append($"MIN_PEP_LENGTH: {MIN_PEP_LENGTH}\n");
            sb.Append($"MAX_PEP_LENGTH: {MAX_PEP_LENGTH}\n");
            sb.Append($"MAX_CHARGE: {MAX_CHARGE}\n");
            sb.Append($"MAX_NEUTRAL_LOSSES: {MAX_NEUTRAL_LOSSES}\n");
            sb.Append($"MAX_NEUTRAL_LOSS_MODS: {MAX_NEUTRAL_LOSS_MODS}\n");
            sb.Append($"MAX_ALLOWED_CHARGE: {MAX_ALLOWED_CHARGE}\n");
            sb.Append($"TOP_N: {TOP_N}\n");
            sb.Append($"TOLERANCE: {TOLERANCE}\n");
            sb.Append($"NORMALIZE: {NORMALIZE}\n");
            sb.Append($"USE_GAUSSIAN: {USE_GAUSSIAN}\n");
            sb.Append("---- SETTINGS END ----");

            return sb.ToString();
        }
    }

    public static class SettingsReader
    {
        public static Settings readSettings(string filename)
        {
            var settings = new Settings();

            using (StreamReader sr = new StreamReader(filename))
            {
                while (!sr.EndOfStream)
                {
                    var line = sr.ReadLine();

                    if (line != null && line.StartsWith("MAX_CLEAVAGES"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MAX_CLEAVAGES = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MIN_PEP_LENGTH"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MIN_PEP_LENGTH = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MAX_PEP_LENGTH"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MAX_PEP_LENGTH = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MAX_CHARGE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MAX_CHARGE = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MAX_NEUTRAL_LOSSES"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MAX_NEUTRAL_LOSSES = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MAX_NEUTRAL_LOSS_MODS"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MAX_NEUTRAL_LOSS_MODS = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MAX_ALLOWED_CHARGE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            settings.MAX_ALLOWED_CHARGE = values[1].Trim();
                        }
                    }

                    if (line != null && line.StartsWith("TOP_N"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.TOP_N = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("TOLERANCE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = float.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.TOLERANCE = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("NORMALIZE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = bool.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.NORMALIZE = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("USE_GAUSSIAN"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = bool.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.USE_GAUSSIAN = value;
                            }
                        }
                    }
                }
            }

            return settings;
        }
    }
}

using System.Text;

namespace CandidateSearch.util
{
    /// <summary>
    /// Settings for digestion, ion calculation and VectorSearch.
    /// </summary>
    public class Settings
    {
        // config digestion
        /// <summary>
        /// Maximum number of missed cleavages allowed during digestion.
        /// </summary>
        public int MAX_CLEAVAGES { get; set; } // default: 2
        /// <summary>
        /// Minimum peptide length.
        /// </summary>
        public int MIN_PEP_LENGTH { get; set; } // default: 5
        /// <summary>
        /// Maximum peptide length.
        /// </summary>
        public int MAX_PEP_LENGTH { get; set; } // default: 30

        // config ion calculation
        /// <summary>
        /// Maximum considered precursor ion charge.
        /// </summary>
        public int MAX_PRECURSOR_CHARGE { get; set; } // default: 4
        /// <summary>
        /// Maximum considered fragment ion charge.
        /// </summary>
        public string MAX_FRAGMENT_CHARGE { get; set; } // default: "+1";
        /// <summary>
        /// Maximum number of considered neutral losses.
        /// </summary>
        public int MAX_NEUTRAL_LOSSES { get; set; } // default: 1
        /// <summary>
        /// Maximum number of considered neutral loss modifications.
        /// </summary>
        public int MAX_NEUTRAL_LOSS_MODS { get; set; } // default: 2
        /// <summary>
        /// Dictionary for fixed modifications that maps amino acids to their possible modification masses.
        /// </summary>
        public Dictionary<string, double> FIXED_MODIFICATIONS { get; set; } // default: empty Dictionary
        /// <summary>
        /// Dictionary for variable modifications that maps amino acids to their possible modification masses.
        /// </summary>
        public Dictionary<string, double> VARIABLE_MODIFICATIONS { get; set; } // default: empty Dictionary

        // config search parameters
        /// <summary>
        /// Whether or not decoy search should be performed.
        /// </summary>
        public bool DECOY_SEARCH { get; set; } // default: true

        // config vector search
        /// <summary>
        /// The top n candidates that should be returned by the VectorSearch.
        /// </summary>
        public int TOP_N { get; set; } // default: 1000
        /// <summary>
        /// The tolerance used for the VectorSearch.
        /// </summary>
        public float TOLERANCE { get; set; } // default: 0.02f
        /// <summary>
        /// Whether or not scores should be normalized by the VectorSearch.
        /// </summary>
        public bool NORMALIZE { get; set; } // default: false
        /// <summary>
        /// Whether or not peaks should be modelled as gaussian distributions by the VectorSearch.
        /// </summary>
        public bool USE_GAUSSIAN { get; set; } // default: true
        /// <summary>
        /// The search approach used by the VectorSearch.
        /// </summary>
        public string MODE { get; set; } // default: CPU_DV

        /// <summary>
        /// Settings constructor to set the specified search parameters.
        /// </summary>
        /// <param name="MaxCleavages">Maximum number of missed cleavages allowed during digestion.</param>
        /// <param name="MinPepLength">Minimum peptide length.</param>
        /// <param name="MaxPepLength">Maximum peptide length.</param>
        /// <param name="MaxPrecursorCharge">Maximum considered precursor ion charge.</param>
        /// <param name="MaxFragmentCharge">Maximum considered fragment ion charge.</param>
        /// <param name="MaxNeutralLosses">Maximum number of considered neutral losses.</param>
        /// <param name="MaxNeutralLossMods">Maximum number of considered neutral loss modifications.</param>
        /// <param name="DecoySearch">Whether or not decoy search should be performed.</param>
        /// <param name="TopN">The top n candidates that should be returned by the VectorSearch.</param>
        /// <param name="Tolerance">The tolerance used for the VectorSearch.</param>
        /// <param name="Normalize">Whether or not scores should be normalized by the VectorSearch.</param>
        /// <param name="UseGaussian">Whether or not peaks should be modelled as gaussian distributions by the VectorSearch.</param>
        /// <param name="Mode">The search approach used by the VectorSearch.</param>
        public Settings(int MaxCleavages = 2, int MinPepLength = 5, int MaxPepLength = 30, 
                        int MaxPrecursorCharge = 4, string MaxFragmentCharge = "+1", int MaxNeutralLosses = 1, int MaxNeutralLossMods = 2,
                        bool DecoySearch = true,
                        int TopN = 1000, float Tolerance = 0.02f, bool Normalize = false, bool UseGaussian = true, string Mode = "CPU_SMi32")
        {
            MAX_CLEAVAGES = MaxCleavages;
            MIN_PEP_LENGTH = MinPepLength;
            MAX_PEP_LENGTH = MaxPepLength;

            MAX_PRECURSOR_CHARGE = MaxPrecursorCharge;
            MAX_FRAGMENT_CHARGE = MaxFragmentCharge;
            MAX_NEUTRAL_LOSSES = MaxNeutralLosses;
            MAX_NEUTRAL_LOSS_MODS = MaxNeutralLossMods;
            FIXED_MODIFICATIONS = new Dictionary<string, double>();
            VARIABLE_MODIFICATIONS = new Dictionary<string, double>();

            DECOY_SEARCH = DecoySearch;

            TOP_N = TopN;
            TOLERANCE = Tolerance;
            NORMALIZE = Normalize;
            USE_GAUSSIAN = UseGaussian;
            MODE = Mode;
        }

        /// <summary>
        /// Add a fixed modification to the fixed modification dictionary.
        /// </summary>
        /// <param name="aminoAcid">The amino acid that will be modified.</param>
        /// <param name="mass">The mass of the modification.</param>
        /// <returns>True if there is no modification for that amino acid yet, false if there already exists a modification for that amino acid.</returns>
        public bool addFixedModification(string aminoAcid, double mass)
        {
            if (FIXED_MODIFICATIONS.ContainsKey(aminoAcid))
                return false;

            FIXED_MODIFICATIONS.Add(aminoAcid, mass);
            return true;
        }

        /// <summary>
        /// Add a variable modification to the variable modification dictionary.
        /// </summary>
        /// <param name="aminoAcid">The amino acid that can be modified.</param>
        /// <param name="mass">The mass of the modification.</param>
        /// <returns>True if there is no modification for that amino acid yet, false if there already exists a modification for that amino acid.</returns>
        public bool addVariableModification(string aminoAcid, double mass)
        {
            if (VARIABLE_MODIFICATIONS.ContainsKey(aminoAcid))
                return false;

            VARIABLE_MODIFICATIONS.Add(aminoAcid, mass);
            return true;
        }

        /// <summary>
        /// Returns a string representation of the modifications.
        /// </summary>
        /// <param name="variable">Whether to process fixed or variable modifications.</param>
        /// <returns>The string representation of the specified modification set.</returns>
        public string modificationsToString(bool variable = false)
        {
            var sb = new StringBuilder();
            if (!variable)
            {
                foreach (var mod in FIXED_MODIFICATIONS)
                {
                    sb.Append($"{mod.Key.ToString()}:{mod.Value.ToString()};");
                }
            }
            else
            {
                foreach (var mod in VARIABLE_MODIFICATIONS)
                {
                    sb.Append($"{mod.Key.ToString()}:{mod.Value.ToString()};");
                }
            }
            return sb.ToString();
        }

        /// <summary>
        /// Returns a string representation of the settings.
        /// </summary>
        /// <returns>The string representation of the settings.</returns>
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("---- SETTINGS: ----\n");
            sb.Append($"MAX_CLEAVAGES: {MAX_CLEAVAGES}\n");
            sb.Append($"MIN_PEP_LENGTH: {MIN_PEP_LENGTH}\n");
            sb.Append($"MAX_PEP_LENGTH: {MAX_PEP_LENGTH}\n");
            sb.Append($"MAX_PRECURSOR_CHARGE: {MAX_PRECURSOR_CHARGE}\n");
            sb.Append($"MAX_FRAGMENT_CHARGE: {MAX_FRAGMENT_CHARGE}\n");
            sb.Append($"MAX_NEUTRAL_LOSSES: {MAX_NEUTRAL_LOSSES}\n");
            sb.Append($"MAX_NEUTRAL_LOSS_MODS: {MAX_NEUTRAL_LOSS_MODS}\n");
            sb.Append($"FIXED_MODIFICATIONS: {modificationsToString(false)}\n");
            sb.Append($"VARIABLE_MODIFICATIONS: {modificationsToString(true)}\n");
            sb.Append($"DECOY_SEARCH: {DECOY_SEARCH}\n");
            sb.Append($"TOP_N: {TOP_N}\n");
            sb.Append($"TOLERANCE: {TOLERANCE}\n");
            sb.Append($"NORMALIZE: {NORMALIZE}\n");
            sb.Append($"USE_GAUSSIAN: {USE_GAUSSIAN}\n");
            sb.Append($"MODE: {MODE}\n");
            sb.Append("---- SETTINGS END ----");

            return sb.ToString();
        }
    }

    /// <summary>
    /// Reader class to read a settings file.
    /// </summary>
    public static class SettingsReader
    {
        /// <summary>
        /// Reads a settings file and returns a settings instance with the adjusted parameters.
        /// </summary>
        /// <param name="filename">The filename of the settings file.</param>
        /// <returns>The settings instance generated from the specified settings file.</returns>
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

                    if (line != null && line.StartsWith("MAX_PRECURSOR_CHARGE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = int.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.MAX_PRECURSOR_CHARGE = value;
                            }
                        }
                    }

                    if (line != null && line.StartsWith("MAX_FRAGMENT_CHARGE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var maxFragmentCharge = values[1].Trim();
                            if (maxFragmentCharge == "+2" ||
                                maxFragmentCharge == "+3" ||
                                maxFragmentCharge == "+4" ||
                                maxFragmentCharge == "Precursor - 1")
                            {
                                settings.MAX_FRAGMENT_CHARGE = maxFragmentCharge;
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

                    if (line != null && line.StartsWith("FIXED_MODIFICATIONS"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var mods = values[1].Split(";");
                            foreach (var mod in mods)
                            {
                                var modProps = mod.Split(":");
                                if (modProps.Length == 2)
                                {
                                    if (modProps[0].Trim().Length == 1)
                                    {
                                        var ok = double.TryParse(modProps[1], 
                                                                 System.Globalization.NumberStyles.AllowDecimalPoint, 
                                                                 System.Globalization.CultureInfo.InvariantCulture, 
                                                                 out var value);
                                        if (ok)
                                        {
                                            if (!settings.FIXED_MODIFICATIONS.ContainsKey(modProps[0]))
                                            {
                                                settings.addFixedModification(modProps[0].Trim(), value);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (line != null && line.StartsWith("VARIABLE_MODIFICATIONS"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var mods = values[1].Split(";");
                            foreach (var mod in mods)
                            {
                                var modProps = mod.Split(":");
                                if (modProps.Length == 2)
                                {
                                    if (modProps[0].Trim().Length == 1)
                                    {
                                        var ok = double.TryParse(modProps[1], 
                                                                 System.Globalization.NumberStyles.AllowDecimalPoint, 
                                                                 System.Globalization.CultureInfo.InvariantCulture, 
                                                                 out var value);
                                        if (ok)
                                        {
                                            if (!settings.VARIABLE_MODIFICATIONS.ContainsKey(modProps[0]))
                                            {
                                                settings.addVariableModification(modProps[0].Trim(), value);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (line != null && line.StartsWith("DECOY_SEARCH"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var ok = bool.TryParse(values[1], out var value);
                            if (ok)
                            {
                                settings.DECOY_SEARCH = value;
                            }
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
                            var ok = float.TryParse(values[1], 
                                                    System.Globalization.NumberStyles.AllowDecimalPoint,
                                                    System.Globalization.CultureInfo.InvariantCulture,
                                                    out var value);
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

                    if (line != null && line.StartsWith("MODE"))
                    {
                        var values = line.Split("=");
                        if (values.Length > 1)
                        {
                            var mode = values[1].Trim();
                            switch (mode)
                            {
                                case "CPU_DV":
                                    settings.MODE = "CPU_DVi32";
                                    break;
                                case "CPU_DVi32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_DVf32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_SV":
                                    settings.MODE = "CPU_SVi32";
                                    break;
                                case "CPU_SVi32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_SVf32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_DM":
                                    settings.MODE = "CPU_DMi32";
                                    break;
                                case "CPU_DMi32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_DMf32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_SM":
                                    settings.MODE = "CPU_SMi32";
                                    break;
                                case "CPU_SMi32":
                                    settings.MODE = mode;
                                    break;
                                case "CPU_SMf32":
                                    settings.MODE = mode;
                                    break;
                                case "GPU_DV":
                                    settings.MODE = "GPU_DVf32";
                                    break;
                                case "GPU_DVf32":
                                    settings.MODE = mode;
                                    break;
                                case "GPU_DM":
                                    settings.MODE = "GPU_DMf32";
                                    break;
                                case "GPU_DMf32":
                                    settings.MODE = mode;
                                    break;
                                case "GPU_SM":
                                    settings.MODE = "GPU_SMf32";
                                    break;
                                case "GPU_SMf32":
                                    settings.MODE = mode;
                                    break;
                                default:
                                    settings.MODE = "CPU_SMi32";
                                    break;
                            }
                        }
                    }
                }
            }

            return settings;
        }
    }
}

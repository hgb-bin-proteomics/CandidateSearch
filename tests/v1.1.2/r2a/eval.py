#!/usr/bin/env python3

import re
import sys
import pandas as pd

def get_candidates(filename, top_n):
    candidates = dict()

    with open(filename, "r", encoding = "utf-8") as f:
        lines = f.readlines()
        f.close()

    for line in lines[1:]:
        values = line.split(";")
        candidates[int(values[0])] = [re.sub(r"[^A-Z]", "", val) for val in values[1].split("],")[:top_n]]

    return candidates

def get_coverage(csms, candidates, top_n, both = False):
    found_candidates = 0

    for csm in csms:
        candidates_for_scan = candidates[csm][:top_n]
        pep1, pep2 = csms[csm]
        if both:
            if pep1 in candidates_for_scan and pep2 in candidates_for_scan:
                found_candidates += 1
        else:
            if pep1 in candidates_for_scan or pep2 in candidates_for_scan:
                found_candidates += 1

    return round((found_candidates/len(csms)) * 100.0, 2)

if __name__ == "__main__":
    maxlynx_csms = pd.read_csv("crosslinkMsms.txt", sep = "\t")
    print(f"Read {maxlynx_csms.shape} CSMs.")
    maxlynx_csms = maxlynx_csms.loc[maxlynx_csms["Raw file"] == f"XLpeplib_Beveridge_QEx-HFX_DSS_R{int(sys.argv[2].strip())}"]
    print(f"Got {maxlynx_csms.shape} CSMs for R{int(sys.argv[2].strip())}.")
    scores = maxlynx_csms["Score"].tolist()
    scores.sort()
    cutoff = 0
    for score in scores[::-1]:
        nr_decoys = maxlynx_csms.loc[(maxlynx_csms["Score"] >= score) & (maxlynx_csms["Decoy"] != "forward")].shape[0]
        nr_targets = maxlynx_csms.loc[(maxlynx_csms["Score"] >= score) & (maxlynx_csms["Decoy"] == "forward")].shape[0]
        if nr_decoys / nr_targets > 0.01:
            cutoff = score
            break
    print(f"Found cutoff for 1% FDR: {cutoff}")
    csms = dict()
    for i, row in maxlynx_csms.iterrows():
        if row["Decoy"].strip() == "forward" and row["Score"] > cutoff:
            csms[int(row["Scan number"])] = [str(row["Sequence1"]).upper(), str(row["Sequence2"]).upper()]
    print(f"Number of CSMs at 1% FDR: {len(csms)}")
    candidates = get_candidates(sys.argv[1], 1000)
    result = pd.DataFrame({"top_n": [1000, 500, 100, 50],
                           "coverage_perc": [get_coverage(csms, candidates, 1000),
                                             get_coverage(csms, candidates, 500),
                                             get_coverage(csms, candidates, 100),
                                             get_coverage(csms, candidates, 50)]})
    result.to_excel(sys.argv[3])
    print("Finished!")

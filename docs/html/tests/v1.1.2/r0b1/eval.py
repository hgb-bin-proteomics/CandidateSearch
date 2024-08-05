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

def get_coverage(psms, candidates, top_n):
    found_candidates = 0

    for psm in psms:
        candidates_for_scan = candidates[psm][:top_n]
        pep = psms[psm]
        if pep in candidates_for_scan:
            found_candidates += 1

    return round((found_candidates/len(psms)) * 100.0, 2)

if __name__ == "__main__":
    amanda_psms = pd.read_excel(sys.argv[2])
    psms = dict()
    for i, row in amanda_psms.iterrows():
        psms[int(row["First Scan"])] = str(row["Annotated Sequence"]).upper()
    candidates = get_candidates(sys.argv[1], 1000)
    result = pd.DataFrame({"top_n": [1000, 500, 100, 50],
                           "coverage_perc": [get_coverage(psms, candidates, 1000),
                                             get_coverage(psms, candidates, 500),
                                             get_coverage(psms, candidates, 100),
                                             get_coverage(psms, candidates, 50)]})
    result.to_excel(sys.argv[3])
    print("Finished!")

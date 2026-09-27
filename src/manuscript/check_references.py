"""
check_references.py - verify every DOI of references.py against Crossref
(the DOI must resolve and the Crossref title must match the cited title).

usage: python src/manuscript/check_references.py
"""
import re
import sys
import time
import unicodedata
from difflib import SequenceMatcher
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from references import REFS  # noqa: E402

DOI = re.compile(r"https://doi\.org/(\S+)$")


def norm(s):
    s = re.sub(r"<[^>]+>", "", s)                   # Crossref titles carry HTML (<sub>, <i>, <scp>)
    s = unicodedata.normalize("NFKD", re.sub(r"[_^]\{([^}]*)\}", r"\1", s))
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    bad = 0
    for key, ref in REFS.items():
        m = DOI.search(ref)
        if not m:
            print(f"  --   {key}: no DOI (checked by hand)")
            continue
        doi = m.group(1)
        after_year = norm(re.split(r"\(\d{4}\) ", ref, maxsplit=1)[1])
        for k in range(4):
            try:
                r = requests.get(f"https://api.crossref.org/works/{doi}", timeout=30,
                                 headers={"User-Agent": "reference-check (mailto:andres.monreal@ues.mx)"})
                break
            except requests.RequestException:
                time.sleep(3)
        if r.status_code != 200:
            if "zenodo" in doi:
                print(f"  ok?  {key}: {doi} (DataCite DOI, not in Crossref)")
                continue
            print(f"  BAD  {key}: {doi} -> HTTP {r.status_code}")
            bad += 1
            continue
        t = (r.json()["message"].get("title") or [""])[0]
        nt = norm(t)
        sim = SequenceMatcher(None, after_year[:len(nt)], nt).ratio()
        if sim <= 0.8 and len(nt) > 60 and nt[:60] in after_year:
            sim = 1.0          # Crossref title = our title + appended footnote (e.g. "Part II: ...")
        flag = "ok " if sim > 0.8 else "BAD"
        bad += flag == "BAD"
        print(f"  {flag}  {key}: {sim:.2f}  {t[:90]}")
        time.sleep(0.2)
    print(f"\n{bad} problem(s)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path

def itd(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0: return None
    rows = list(csv.DictReader(p.open(), delimiter="\t"))
    return rows[0] if rows else None

def pizzly_call(path, a, b):
    data = json.load(open(path))
    for record in data.get("genes", []):
        observed = {record["geneA"]["name"], record["geneB"]["name"]}
        if observed == {a, b}:
            return record
    return None

ap = argparse.ArgumentParser()
ap.add_argument("--mv-itd", required=True); ap.add_argument("--kasumi-itd", required=True)
ap.add_argument("--kasumi-fusion", required=True); ap.add_argument("--nb4-fusion", required=True)
ap.add_argument("--kasumi-pizzly", required=True); ap.add_argument("--nb4-pizzly", required=True)
ap.add_argument("--outdir", required=True); a=ap.parse_args()
out=Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
mv=itd(a.mv_itd); ka=itd(a.kasumi_itd)
kp=pizzly_call(a.kasumi_pizzly,"RUNX1","RUNX1T1")
np=pizzly_call(a.nb4_pizzly,"PML","RARA")
rows=[
 ["MV4-11","FLT3-ITD","positive","positive" if mv else "negative","PASS" if mv else "FAIL", mv.get("length","") if mv else "",mv.get("vaf","") if mv else "",mv.get("counts_unique_total","") if mv else ""],
 ["KASUMI-1","FLT3-ITD","negative","positive" if ka else "negative","PASS" if not ka else "FAIL",ka.get("length","") if ka else "",ka.get("vaf","") if ka else "",ka.get("counts_unique_total","") if ka else ""],
 ["KASUMI-1","RUNX1--RUNX1T1","positive","positive" if kp else "negative","PASS" if kp else "FAIL","","",f"pizzly split={kp.get('splitcount',0) if kp else 0};pair={kp.get('paircount',0) if kp else 0}"],
 ["NB-4","PML--RARA","positive","positive" if np else "negative","PASS" if np else "FAIL","","",f"pizzly split={np.get('splitcount',0) if np else 0};pair={np.get('paircount',0) if np else 0}"],]
with (out/"known_positive_recovery.tsv").open("w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["sample","event","expected","observed","result","length_bp","vaf_percent","evidence"]); w.writerows(rows)
with (out/"summary.txt").open("w") as fh:
    fh.write("Track 2 known-positive recovery (not a sensitivity/specificity validation)\n")
    for r in rows: fh.write(f"{r[0]} {r[1]}: expected {r[2]}, observed {r[3]} [{r[4]}]; evidence {r[7] or 'see fields'}\n")
    if kp and np:
        fh.write("Fusion cells PASS: pizzly recovered both preregistered expected gene pairs.\n")
    else:
        fh.write("At least one preregistered expected fusion was not recovered by pizzly.\n")

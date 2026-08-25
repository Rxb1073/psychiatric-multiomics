import csv, gzip, os, re
from collections import defaultdict
PICR = "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/picrust2"
OUT  = os.path.join(PICR, "pathway_DA")
DATASETS = ["PRJNA591924", "PRJNA776170"]
DEPLETED = ["PWY490.3","PWY.6609","OANTIGEN.PWY","PWY4FS.7","PWY4FS.8","PWY0.1297",
            "TRNA.CHARGING.PWY","P161.PWY","UDPNAGSYN.PWY","GLYCOGENSYNTH.PWY",
            "PWY.5100","PWY.6126","ANAGLYCOLYSIS.PWY","PHOSLIPSYN.PWY"]
BUTYRATE = ["CENTFERM.PWY","PWY.5676","P162.PWY","PWY.6590","PWY.5022"]
def variants(pid): return {pid, pid.replace(".","-")}
TARGET_MATCH = {}
for p in DEPLETED:
    for v in variants(p): TARGET_MATCH[v] = (p, "depleted")
for p in BUTYRATE:
    for v in variants(p): TARGET_MATCH[v] = (p, "butyrate")
SCFA_GENERA = {"Faecalibacterium","Fusicatenibacter","Barnesiella"}
SCFA_PARTIAL = ["Lachnospiraceae"]
def load_tax(d):
    m = {}
    with open(os.path.join(PICR, d, "taxonomy_export", "taxonomy.tsv")) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 2: continue
            mm = re.search(r"g__([^;]*)", f[1])
            m[f[0]] = (mm.group(1).strip() if mm and mm.group(1).strip() else "unclassified")
    return m
def process(d, gc):
    tax = load_tax(d)
    p = os.path.join(PICR, d, "picrust2_output", "pathways_out", "path_abun_contrib.tsv.gz")
    with gzip.open(p, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            hit = TARGET_MATCH.get(row["function"])
            if not hit: continue
            pid, kind = hit
            genus = tax.get(row["taxon"], "unclassified")
            try: val = float(row["taxon_function_abun"])
            except ValueError: continue
            gc[(kind, genus)] += val
def is_scfa(g):
    return g in SCFA_GENERA or any(part in g for part in SCFA_PARTIAL)
gc = defaultdict(float)
for d in DATASETS:
    print(f"processing {d} ...", flush=True); process(d, gc)
def report(kind, label):
    rows = [(g, v) for (k,g), v in gc.items() if k==kind]
    total = sum(v for _,v in rows); rows.sort(key=lambda x: -x[1])
    print(f"\n===== {label}: top genus contributors =====")
    print(f"{'genus':35s} {'contribution':>14s} {'share%':>8s} {'SCFA?':>6s}")
    for g,v in rows[:15]:
        sh = 100*v/total if total else 0
        print(f"{g:35s} {v:14.1f} {sh:8.2f} {'YES' if is_scfa(g) else '':>6s}")
    scfa_total = sum(v for g,v in rows if is_scfa(g))
    print(f"\n  --> SCFA-producer genera = {100*scfa_total/total:.1f}% of {label} contribution" if total else "  none")
report("depleted", "DEPLETED-in-MDD pathways")
report("butyrate", "BUTYRATE/SCFA pathways")
with open(os.path.join(OUT, "taxonomy_function_bridge.tsv"), "w") as fh:
    fh.write("kind\tgenus\tcontribution\tshare_pct\tis_scfa_producer\n")
    for kind in ("depleted","butyrate"):
        rows = [(g,v) for (k,g),v in gc.items() if k==kind]
        total = sum(v for _,v in rows); rows.sort(key=lambda x:-x[1])
        for g,v in rows:
            fh.write(f"{kind}\t{g}\t{v:.1f}\t{100*v/total if total else 0:.2f}\t{int(is_scfa(g))}\n")
print("\nWROTE taxonomy_function_bridge.tsv")

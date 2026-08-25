#!/usr/bin/env python3
# ASD taxonomy-function bridge — replicates SCZ bridge logic across ASD's 4 datasets
import csv, gzip, os, re
from collections import defaultdict

BASE="/rds/projects/e/elhamsak-dissertation-autism/Microbiome"
OUT=os.path.join(BASE,"pooled")
DATASETS={
 "PRJNA453621":  ("PRJNA453621/picrust2_out/pathways_out/path_abun_contrib.tsv.gz","PRJNA453621/tax_export/taxonomy.tsv"),
 "PRJNA529598":  ("PRJNA529598/picrust2_out/pathways_out/path_abun_contrib.tsv.gz","PRJNA529598/tax_export/taxonomy.tsv"),
 "PRJNA769228":  ("PRJNA769228/picrust2_out/pathways_out/path_abun_contrib.tsv.gz","PRJNA769228/tax_export/taxonomy.tsv"),
 "PRJNA1022510": ("PRJNA1022510/picrust2_out/pathways_out/path_abun_contrib.tsv.gz","PRJNA1022510/tax_export/taxonomy.tsv"),
}
PATHDA=os.path.join(OUT,"stage6_pathway_results.csv")

dep=set()
with open(PATHDA) as f:
    r=csv.DictReader(f)
    for row in r:
        try:
            if float(row['coef'])<0 and float(row['qval'])<0.05:
                pid=row['feature']; dep|={pid, pid.replace('.','-'), pid.replace('-','.')}
        except: pass
print("depleted pathways (coef<0,q<0.05):", len(set(p for p in dep if '-' not in p)))

BUT_raw=["CENTFERM.PWY","FERMENTATION.PWY","PWY.5676","P162.PWY","PWY.6590","PWY.5022","PWY.5088"]
BUT=set()
for p in BUT_raw: BUT|={p, p.replace(".","-")}

SCFA={"Faecalibacterium","Fusicatenibacter","Anaerostipes","Blautia","Coprococcus",
      "Roseburia","Agathobacter","Subdoligranulum","Butyricicoccus","Marvinbryantia"}
SCFA_PART=["Lachnospiraceae"]
def is_scfa(g): return g in SCFA or any(p in g for p in SCFA_PART)

def load_tax(path):
    m={}
    with open(path) as f:
        r=csv.reader(f,delimiter='\t'); next(r)
        for row in r:
            if len(row)<2: continue
            mm=re.search(r"g__([^;]*)", row[1])
            m[row[0]]=(mm.group(1).strip() if mm and mm.group(1).strip() else "unclassified")
    return m

dep_c=defaultdict(float); but_c=defaultdict(float)
for d,(contrib,taxf) in DATASETS.items():
    cp=os.path.join(BASE,contrib); tp=os.path.join(BASE,taxf)
    if not os.path.exists(cp): print(f"  SKIP {d}: no contrib file"); continue
    tax=load_tax(tp)
    print(f"processing {d} ...", flush=True)
    with gzip.open(cp,'rt') as f:
        r=csv.DictReader(f,delimiter='\t')
        for row in r:
            fn=row['function']; g=tax.get(row['taxon'],'unclassified')
            try: v=float(row['taxon_function_abun'])
            except: continue
            if fn in dep: dep_c[g]+=v
            if fn in BUT: but_c[g]+=v

def report(gc,label):
    rows=sorted(gc.items(),key=lambda x:-x[1]); tot=sum(v for _,v in rows)
    print(f"\n===== {label}: top genus contributors =====")
    print(f"{'genus':35s} {'share%':>8s} {'SCFA?':>6s}")
    for g,v in rows[:15]:
        print(f"{g:35s} {100*v/tot if tot else 0:8.2f} {'YES' if is_scfa(g) else '':>6s}")
    s=sum(v for g,v in rows if is_scfa(g))
    print(f"  --> SCFA producers = {100*s/tot if tot else 0:.1f}% of {label}")
    return rows,tot

dr,dt=report(dep_c,f"DEPLETED-in-ASD pathways (n={len(set(p for p in dep if '-' not in p))})")
br,bt=report(but_c,"BUTYRATE/SCFA pathways")

with open(os.path.join(OUT,"taxonomy_function_bridge_asd.tsv"),"w") as w:
    w.write("kind\tgenus\tcontribution\tshare_pct\tis_scfa\n")
    for kind,(rows,tot) in [("depleted",(dr,dt)),("butyrate",(br,bt))]:
        for g,v in rows:
            w.write(f"{kind}\t{g}\t{v:.1f}\t{100*v/tot if tot else 0:.2f}\t{int(is_scfa(g))}\n")
print("\nWROTE taxonomy_function_bridge_asd.tsv")

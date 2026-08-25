#!/usr/bin/env python3
import csv, gzip, os
from collections import defaultdict
base="/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407"
contrib=os.path.join(base,"picrust2/picrust2_output/pathways_out/path_abun_contrib.tsv.gz")
pathDA=os.path.join(base,"picrust2/pathway_DA/pathway_SZ_vs_HC.tsv")
taxf=os.path.join(base,"qiime2_v4/_tax/taxonomy.tsv")
OUT=os.path.join(base,"picrust2/pathway_DA")

# depleted pathways (coef<0, q<0.05) pulled dynamically
dep=set()
with open(pathDA) as f:
    r=csv.DictReader(f,delimiter='\t')
    for row in r:
        try:
            if float(row['coef'])<0 and float(row['qval'])<0.05:
                pid=row['feature']; dep.add(pid); dep.add(pid.replace('.','-'))
        except: pass
# butyrate/SCFA-specific pathway set (both . and - variants)
BUT_raw=["CENTFERM.PWY","FERMENTATION.PWY","PWY.5676","P162.PWY","PWY.6590","PWY.5022","PWY.5088"]
BUT=set()
for p in BUT_raw: BUT|={p, p.replace(".","-")}
print("depleted pathways:", len(dep))

# ASV -> genus
asv2gen={}
with open(taxf) as f:
    r=csv.reader(f,delimiter='\t'); next(r)
    for row in r:
        import re
        mm=re.search(r"g__([^;]*)", row[1])
        asv2gen[row[0]] = (mm.group(1).strip() if mm and mm.group(1).strip() else "unclassified")

SCFA={"Faecalibacterium","Fusicatenibacter","Anaerostipes","Blautia","Coprococcus",
      "Roseburia","Agathobacter","Subdoligranulum","Butyricicoccus","Marvinbryantia"}
SCFA_PART=["Lachnospiraceae"]
def is_scfa(g): return g in SCFA or any(p in g for p in SCFA_PART)

dep_c=defaultdict(float); but_c=defaultdict(float)
with gzip.open(contrib,'rt') as f:
    r=csv.DictReader(f,delimiter='\t')
    for row in r:
        fn=row['function']; g=asv2gen.get(row['taxon'],'unclassified')
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
dr,dt=report(dep_c,f"DEPLETED-in-SZ pathways (n={len(dep)})")
br,bt=report(but_c,"BUTYRATE/SCFA pathways")

with open(os.path.join(OUT,"taxonomy_function_bridge_scz.tsv"),"w") as w:
    w.write("kind\tgenus\tcontribution\tshare_pct\tis_scfa\n")
    for kind,(rows,tot) in [("depleted",(dr,dt)),("butyrate",(br,bt))]:
        for g,v in rows:
            w.write(f"{kind}\t{g}\t{v:.1f}\t{100*v/tot if tot else 0:.2f}\t{int(is_scfa(g))}\n")
print("\nWROTE taxonomy_function_bridge_scz.tsv")

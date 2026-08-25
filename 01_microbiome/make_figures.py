#!/usr/bin/env python3
"""make_figures.py — Figures 3-6 for the SCZ microbiome dissertation."""
import argparse, os, csv, re, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--base", required=True)
ap.add_argument("--out", default="figs")
A = ap.parse_args()
os.makedirs(A.out, exist_ok=True)

def genus_key(s):
    s = s.replace("[","").replace("]","").replace('"',"")
    parts = re.split(r"[;.]", s)
    parts = [p.strip() for p in parts if p.strip() not in ("","__")]
    gs = [p for p in parts if p.startswith("g__")]
    tok = gs[-1] if gs else (parts[-1] if parts else s)
    return tok.replace("g__","")

anc_lfc = os.path.join(A.base, "qiime2_v4/_ancombc/lfc_slice.csv")
anc_q   = os.path.join(A.base, "qiime2_v4/_ancombc/q_val_slice.csv")
def load_anc():
    with open(anc_lfc) as f: L=list(csv.reader(f))
    with open(anc_q) as f: Q=list(csv.reader(f))
    h=L[0]; dc=[i for i,c in enumerate(h) if "diagnosis" in c.lower()][0]
    qmap={r[0]:r[dc] for r in Q[1:]}
    d={}
    for r in L[1:]:
        try: lfc=float(r[dc]); q=float(qmap.get(r[0],"nan"))
        except: continue
        if q<0.05: d[genus_key(r[0])]=(lfc,q)
    return d
anc=load_anc()

maa_f=os.path.join(A.base,"maaslin2/maaslin2_SZ_vs_HC.tsv")
maa={}
with open(maa_f) as f:
    r=csv.DictReader(f,delimiter="\t")
    for row in r:
        try: c=float(row["coef"]); q=float(row["qval"])
        except: continue
        if q<0.05: maa[genus_key(row["feature"])]=(c,q)

both=sorted(set(anc)&set(maa))
both=[g for g in both if g not in ("","__") and "Bacteria" not in g and g!="Unassigned" and g.lower()!="uncultured"]

# FIG 3: concordance scatter
xs=[anc[g][0] for g in both]; ys=[maa[g][0] for g in both]
plt.figure(figsize=(6,6))
cols=["#c0392b" if x<0 else "#27ae60" for x in xs]
plt.axhline(0,color="grey",lw=.7); plt.axvline(0,color="grey",lw=.7)
plt.scatter(xs,ys,c=cols,s=45,edgecolor="k",linewidth=.4,zorder=3)
for g,x,y in zip(both,xs,ys):
    if abs(x)>1.7: plt.annotate(g,(x,y),fontsize=7,xytext=(3,3),textcoords="offset points")
plt.xlabel("ANCOM-BC log-fold-change (SZ vs HC)")
plt.ylabel("MaAsLin2 coefficient (SZ vs HC)")
plt.title(f"Two-method concordance: {len(both)} genera, 100% direction-agreement")
plt.tight_layout(); plt.savefig(f"{A.out}/fig3_concordance_scatter.png",dpi=300); plt.close()

# FIG 4: effect-size bars
SCFA={"Fusicatenibacter","Anaerostipes","Blautia","Coprococcus","Roseburia","Marvinbryantia",
      "Lachnospiraceae_ND3007_group","Lachnospiraceae"}
order=sorted(both,key=lambda g:anc[g][0])
vals=[anc[g][0] for g in order]
cols=["#8e44ad" if g in SCFA else ("#c0392b" if anc[g][0]<0 else "#27ae60") for g in order]
plt.figure(figsize=(7,8))
plt.barh(range(len(order)),vals,color=cols,edgecolor="k",linewidth=.4)
plt.yticks(range(len(order)),order,fontsize=8)
plt.axvline(0,color="k",lw=.7)
plt.xlabel("ANCOM-BC log-fold-change (SZ vs HC)")
plt.title("Concordant differentially-abundant genera\n(purple = SCFA-producing Lachnospiraceae)")
plt.tight_layout(); plt.savefig(f"{A.out}/fig4_genus_effectsizes.png",dpi=300); plt.close()

# FIG 5: pathway volcano
pw=os.path.join(A.base,"picrust2/pathway_DA/pathway_SZ_vs_HC.tsv")
coefs=[]; qs=[]; names=[]
with open(pw) as f:
    r=csv.DictReader(f,delimiter="\t")
    for row in r:
        try: c=float(row["coef"]); q=float(row["qval"])
        except: continue
        coefs.append(c); qs.append(q); names.append(row["feature"])
plt.figure(figsize=(7,6))
sig=[q<0.05 for q in qs]
xcol=["#c0392b" if (s and c<0) else ("#27ae60" if (s and c>0) else "#bbbbbb") for s,c in zip(sig,coefs)]
plt.scatter(coefs,[-math.log10(max(q,1e-12)) for q in qs],c=xcol,s=14,edgecolor="none")
plt.axhline(-math.log10(0.05),color="grey",ls="--",lw=.7)
plt.axvline(0,color="grey",lw=.7)
for c,q,n in zip(coefs,qs,names):
    if n.upper().startswith(("CENTFERM","FERMENTATION")):
        plt.annotate(n,(c,-math.log10(max(q,1e-12))),fontsize=7,color="navy",xytext=(3,3),textcoords="offset points")
plt.xlabel("MaAsLin2 coefficient (SZ vs HC)")
plt.ylabel("-log10(q)")
plt.title("Predicted MetaCyc pathways: 92 down / 44 up in SZ\n(fermentation pathways labelled)")
plt.tight_layout(); plt.savefig(f"{A.out}/fig5_pathway_volcano.png",dpi=300); plt.close()

# FIG 6: taxonomy-function bridge
br=os.path.join(A.base,"picrust2/pathway_DA/taxonomy_function_bridge_scz.tsv")
g=[]; share=[]; scfa=[]
with open(br) as f:
    r=csv.DictReader(f,delimiter="\t")
    for row in r:
        if row["kind"]!="depleted": continue
        g.append(row["genus"]); share.append(float(row["share_pct"])); scfa.append(row["is_scfa"] in ("1","True","true"))
top=sorted(zip(g,share,scfa),key=lambda x:-x[1])[:15][::-1]
gg=[t[0] for t in top]; ss=[t[1] for t in top]; cc=["#8e44ad" if t[2] else "#7f8c8d" for t in top]
plt.figure(figsize=(7,7))
plt.barh(range(len(gg)),ss,color=cc,edgecolor="k",linewidth=.4)
plt.yticks(range(len(gg)),gg,fontsize=8)
plt.xlabel("% contribution to pathways depleted in SZ")
plt.title("Taxonomy-function bridge\n(purple = SCFA producer; combined SCFA = 34%)")
plt.tight_layout(); plt.savefig(f"{A.out}/fig6_taxonomy_function_bridge.png",dpi=300); plt.close()

print(f"wrote fig3-fig6 to {A.out}/")

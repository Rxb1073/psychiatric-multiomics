import pandas as pd, gzip
studies={"PRJNA453621":"../PRJNA453621","PRJNA529598":"../PRJNA529598",
         "PRJNA769228":"../PRJNA769228","PRJNA1022510":"../PRJNA1022510"}
tabs=[]
for st,base in studies.items():
    f=f"{base}/picrust2_out/pathways_out/path_abun_unstrat.tsv.gz"
    df=pd.read_csv(f,sep="\t",index_col=0)
    df.columns=[c.strip() for c in df.columns]
    tabs.append(df)
    print(f"{st}: {df.shape[1]} samples, {df.shape[0]} pathways")

# outer-join on pathway ID (union of pathways, samples side by side), fill missing 0
merged=pd.concat(tabs,axis=1).fillna(0)
# collapse duplicate sample columns if any (shouldn't be)
merged=merged.groupby(level=0,axis=1).sum()
print("MERGED:",merged.shape[1],"samples,",merged.shape[0],"pathways")

# align to pooled metadata (only harmonized samples)
meta=pd.read_csv("metadata_pooled.tsv",sep="\t")
meta=meta[meta['sample-id']!="#q2:types"].set_index('sample-id')
common=[s for s in merged.columns if s in meta.index]
print("samples matching metadata:",len(common),"of",merged.shape[1])
merged=merged[common]
merged.to_csv("pooled_pathways.tsv",sep="\t")
meta.loc[common].to_csv("pooled_pathways_meta.tsv",sep="\t")
print("wrote pooled_pathways.tsv + pooled_pathways_meta.tsv")

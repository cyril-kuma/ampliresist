#!/usr/bin/env python3
"""Roadmap-driven publication figures for DRAG1 molecular xenomonitoring.

The five figures deliberately follow the narrative in
``Bioinformatics Visualization Roadmap Development.md``.  Metadata are read
only from the cleaned CSV supplied with ``--metadata``.  Genomic calls and
callability remain products of the pipeline and are never inferred from
missing VCF records.
"""
from __future__ import annotations

import argparse, json, math, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Rectangle
from matplotlib.lines import Line2D

OI = {"orange":"#E69F00","sky":"#56B4E9","green":"#009E73","yellow":"#F0E442",
      "blue":"#0072B2","vermillion":"#D55E00","purple":"#CC79A7","black":"#111111"}
ZONE_ORDER = ["Coastal_Savannah","Forest","Northern_Savannah"]
ZONE_LABEL = {"Coastal_Savannah":"Coastal Savannah","Forest":"Forest","Northern_Savannah":"Northern Savannah"}
SPECIES_ORDER = ["An_arabiensis","An_coluzzii","An_gambiae_s.s"]
SPECIES_LABEL = {"An_arabiensis":r"$An.$ arabiensis","An_coluzzii":r"$An.$ coluzzii","An_gambiae_s.s":r"$An.$ gambiae s.s."}
MARKERS = {
    "dhfr_152_AT":r"$pfdhfr$ N51I", "dhfr_175_TC":r"$pfdhfr$ C59R",
    "dhfr_323_GA":r"$pfdhfr$ S108N", "dhps_1306_TG":r"$pfdhps$ A437G",
    "dhps_1310_GC":r"$pfdhps$ K540E", "dhps_1837_GT":r"$pfdhps$ A581G",
    "mdr1_551_AT":r"$pfmdr1$ Y184F", "mdr1_256_AT":r"$pfmdr1$ N86Y",
    "crt_227_AC":r"$pfcrt$ K76T"
}

def style():
    mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"axes.titlesize":10,
        "axes.labelsize":9,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5,
        "axes.linewidth":0.7,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
        "savefig.dpi":400,"savefig.bbox":"tight","axes.spines.top":False,"axes.spines.right":False})

def panel(ax, letter, title):
    ax.text(-.10,1.015,letter,transform=ax.transAxes,fontsize=12,fontweight="bold",va="bottom")
    ax.set_title(title,loc="left",fontweight="bold",pad=5)

def export(fig, outdir, stem):
    for ext in ("png","svg","pdf"):
        fig.savefig(outdir/f"{stem}.{ext}",dpi=400 if ext=="png" else None,facecolor="white")
    plt.close(fig)

def wilson(k,n,z=1.959963984540054):
    if not n: return np.nan,np.nan
    p=k/n; den=1+z*z/n; cen=(p+z*z/(2*n))/den; half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,cen-half),min(1,cen+half)

def load_geo(path):
    if not path or not Path(path).exists(): return []
    obj=json.loads(Path(path).read_text())
    rings=[]
    def walk(x):
        if isinstance(x,list) and x and isinstance(x[0],list) and len(x[0])>=2 and all(isinstance(v,(int,float)) for v in x[0][:2]): rings.append(x)
        elif isinstance(x,list):
            for y in x: walk(y)
    for f in obj.get("features",[]): walk(f.get("geometry",{}).get("coordinates",[]))
    return rings

def logistic_ridge(df, outcome):
    """Site-adjusted penalised logistic regression using deterministic IRLS.

    Collection-site indicators receive a weak ridge penalty, approximating
    partial pooling of site intercepts. Fixed-effect estimates are unpenalised.
    Wald intervals are explicitly labelled approximate in figure/caption.
    """
    use=df[[outcome,"bioclimatic_zone","sibling_species","host_feeding_type","pf_ct","collection_site"]].dropna().copy()
    use=use[use.host_feeding_type.isin(["single_host","mixed_host"])]
    refs={"bioclimatic_zone":"Forest","sibling_species":"An_coluzzii","host_feeding_type":"single_host","collection_site":"Prang"}
    cols=[np.ones(len(use))]; names=["Intercept"]
    for var,levels in [("bioclimatic_zone",ZONE_ORDER),("sibling_species",SPECIES_ORDER),("host_feeding_type",["single_host","mixed_host"])]:
        for level in levels:
            if level==refs[var] or level not in set(use[var]): continue
            cols.append((use[var].to_numpy()==level).astype(float)); names.append(f"{var}:{level}")
    ct=use.pf_ct.astype(float).to_numpy(); ct_sd=ct.std() or 1
    cols.append((ct-ct.mean())/ct_sd); names.append("pf_ct:per_SD")
    fixed_n=len(cols)
    for level in sorted(set(use.collection_site)-{refs["collection_site"]}):
        cols.append((use.collection_site.to_numpy()==level).astype(float)); names.append(f"site:{level}")
    X=np.column_stack(cols); y=use[outcome].astype(float).to_numpy(); beta=np.zeros(X.shape[1]); penalty=np.zeros(X.shape[1]); penalty[fixed_n:]=1.0
    for _ in range(200):
        eta=np.clip(X@beta,-25,25); mu=1/(1+np.exp(-eta)); w=np.maximum(mu*(1-mu),1e-7)
        H=X.T@(w[:,None]*X)+np.diag(penalty); g=X.T@(y-mu)-penalty*beta
        step=np.linalg.pinv(H)@g; beta+=step
        if np.max(np.abs(step))<1e-8: break
    cov=np.linalg.pinv(H); se=np.sqrt(np.maximum(np.diag(cov),0))
    rows=[]
    labels={"bioclimatic_zone:Coastal_Savannah":"Zone: Coastal Savannah vs Forest",
      "bioclimatic_zone:Northern_Savannah":"Zone: Northern Savannah vs Forest",
      "sibling_species:An_arabiensis":r"Species: $An.$ arabiensis vs $An.$ coluzzii",
      "sibling_species:An_gambiae_s.s":r"Species: $An.$ gambiae s.s. vs $An.$ coluzzii",
      "host_feeding_type:mixed_host":"Feeding: mixed vs single host","pf_ct:per_SD":f"Pf Ct: per {ct_sd:.1f}-cycle SD"}
    for i,nm in enumerate(names[:fixed_n]):
        if nm=="Intercept": continue
        rows.append({"term":nm,"label":labels.get(nm,nm),"log_or":beta[i],"se":se[i],"odds_ratio":math.exp(np.clip(beta[i],-20,20)),
                     "conf_low":math.exp(np.clip(beta[i]-1.96*se[i],-20,20)),"conf_high":math.exp(np.clip(beta[i]+1.96*se[i],-20,20)),"n":len(use),"events":int(y.sum())})
    return pd.DataFrame(rows), use

def figure1(d, adm0, adm1, out):
    # Site-level DHFR haplotypes communicate geography and clinically meaningful composition.
    haps=[x for x in ["NCSI","NRNI","IRNI","ICNI"] if x in set(d.dhfr_haplotype.dropna())]
    hc={h:c for h,c in zip(haps,["#BDBDBD",OI["sky"],OI["orange"],OI["purple"]])}
    q=d.dropna(subset=["dhfr_haplotype","collection_site","latitude","longitude"])
    ct=q.groupby(["collection_site","latitude","longitude","bioclimatic_zone","dhfr_haplotype"]).size().unstack(fill_value=0).reindex(columns=haps,fill_value=0).reset_index()
    ct["n_callable"]=ct[haps].sum(axis=1); ct.to_csv(out/"figure1_geographic_source.csv",index=False)
    fig,ax=plt.subplots(figsize=(7.2,6.2)); panel(ax,"A","Geographic composition of pfdhfr haplotypes")
    for ring in adm0: ax.fill([p[0] for p in ring],[p[1] for p in ring],fc="#FAFAFA",ec="#555",lw=.8,zorder=0)
    for ring in adm1: ax.plot([p[0] for p in ring],[p[1] for p in ring],c="#D0D0D0",lw=.45,zorder=.5)
    nudges={"Aplaku":(-.18,-.05),"Teshie":(.18,.04),"Obuasi":(-.12,.04),"Bekwai":(.12,.08),"Chiraa":(-.08,.08),"Sunyani":(.10,-.04)}
    for _,r in ct.sort_values("n_callable").iterrows():
        x,y=float(r.longitude),float(r.latitude); dx,dy=nudges.get(r.collection_site,(0,0)); x2,y2=x+dx,y+dy; n=int(r.n_callable); rad=.035+.018*math.sqrt(n)
        if dx or dy: ax.plot([x,x2],[y,y2],c="#777",lw=.5,zorder=2)
        start=90
        for h in haps:
            extent=360*int(r[h])/n
            ax.add_patch(Wedge((x2,y2),rad,start,start+extent,facecolor=hc[h],edgecolor="white",lw=.5,zorder=3)); start+=extent
        ax.add_patch(Circle((x2,y2),rad,fill=False,ec="#333",lw=.5,zorder=4))
        ax.annotate(f"{r.collection_site}\n$n$={n}",(x2,y2),xytext=(rad*1.15,0),textcoords="offset points",fontsize=6.5,va="center",fontweight="bold")
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)"); ax.set_aspect("equal"); ax.grid(False)
    ax.legend([Rectangle((0,0),1,1,fc=hc[h]) for h in haps],haps,title=r"$pfdhfr$ haplotype",frameon=False,loc="upper left",bbox_to_anchor=(1.01,.98))
    ax.text(.01,.01,"Pie area ∝ DHFR-callable specimens; slices show within-site composition.\nSmall denominators are descriptive, not population prevalence.",transform=ax.transAxes,fontsize=7,color="#444")
    fig.suptitle("Figure 1 | Where are pyrimethamine-resistance haplotypes circulating?",x=.08,ha="left",fontsize=13,fontweight="bold")
    fig.subplots_adjust(left=.12,right=.78,bottom=.11,top=.86); export(fig,out,"figure1_geographic_resistance")

def figure2(d,out):
    keys=["dhfr_152_AT","dhfr_175_TC","dhfr_323_GA","dhps_1306_TG","dhps_1310_GC"]
    q=d.dropna(subset=keys).copy(); bits=q[keys].astype(int).astype(str).agg("".join,axis=1); counts=bits.value_counts().head(15)
    src=[]
    for pattern,n in counts.items():
        row={"pattern":pattern,"n":n}; row.update({k:int(v) for k,v in zip(keys,pattern)}); src.append(row)
    src=pd.DataFrame(src); src.to_csv(out/"figure2_upset_source.csv",index=False)
    fig=plt.figure(figsize=(7.2,5.6)); gs=fig.add_gridspec(2,1,height_ratios=[2.2,1],hspace=.04); ax=fig.add_subplot(gs[0]); mx=fig.add_subplot(gs[1],sharex=ax)
    x=np.arange(len(src)); ax.bar(x,src.n,color=OI["blue"],width=.75); ax.set_ylabel("Specimens"); ax.set_xticks([]); panel(ax,"A","Most frequent clinically relevant mutation intersections")
    for i,n in enumerate(src.n): ax.text(i,n+max(src.n)*.015,str(n),ha="center",va="bottom",fontsize=7,fontweight="bold")
    for row,k in enumerate(keys):
        yy=len(keys)-1-row; mx.text(-.75,yy,MARKERS[k],ha="right",va="center",fontsize=8)
        mx.scatter(x,np.full(len(x),yy),s=15,c="#D8D8D8",zorder=1)
        on=np.where(src[k].to_numpy()==1)[0]; mx.scatter(on,np.full(len(on),yy),s=24,c=OI["black"],zorder=2)
        for i in range(len(src)):
            ys=np.where(src.loc[i,keys].to_numpy()==1)[0]
            if len(ys)>1: mx.plot([i,i],[len(keys)-1-ys.max(),len(keys)-1-ys.min()],c=OI["black"],lw=1.2,zorder=1)
    mx.set_ylim(-.7,len(keys)-.3); mx.set_yticks([]); mx.set_xticks(x); mx.set_xticklabels([f"I{i+1}" for i in x],fontsize=7); mx.set_xlabel("Mutation intersection (complete calls only)"); mx.spines[["left","right","top"]].set_visible(False)
    fig.suptitle("Figure 2 | Multigenic architecture of SP resistance",x=.08,ha="left",fontsize=13,fontweight="bold")
    fig.text(.08,.01,f"Complete-case denominator: $n$={len(q)}. Grey = absent; black = present. Missing calls are excluded, never treated as wild type.",fontsize=7,color="#444")
    fig.tight_layout(rect=(.12,.04,1,.94)); export(fig,out,"figure2_multigenic_upset")

def figure3(d,out):
    keys=[k for k in MARKERS if k in d and pd.to_numeric(d[k],errors="coerce").sum(skipna=True)>0]
    rows=[]
    for facet,var,levels in [("Bioclimatic zone","bioclimatic_zone",ZONE_ORDER),("Vector species","sibling_species",SPECIES_ORDER)]:
        for level in levels:
            for k in keys:
                z=pd.to_numeric(d.loc[d[var]==level,k],errors="coerce").dropna(); n=len(z); mut=int((z>0).sum()); lo,hi=wilson(mut,n)
                rows.append({"facet":facet,"group":level,"marker":k,"mutant":mut,"callable":n,"frequency":mut/n if n else np.nan,"low":lo,"high":hi})
    s=pd.DataFrame(rows); s.to_csv(out/"figure3_abacus_source.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.2,4.8),sharey=True,gridspec_kw={"wspace":.12})
    for a,(facet,sub),letter in zip(axs,s.groupby("facet",sort=False),"AB"):
        groups=list(dict.fromkeys(sub.group)); ypos={k:i for i,k in enumerate(keys[::-1])}
        offsets=np.linspace(-.22,.22,len(groups)); colors=[OI["sky"],OI["green"],OI["orange"]]
        for off,g,c in zip(offsets,groups,colors):
            z=sub[sub.group==g]
            yy=np.array([ypos[k] for k in z.marker])+off; xx=100*z.frequency.to_numpy();
            a.hlines(yy,100*z.low,100*z.high,color=c,lw=1); a.scatter(xx,yy,s=18+2*np.sqrt(z.callable),c=c,edgecolor="white",lw=.35,label=ZONE_LABEL.get(g,SPECIES_LABEL.get(g,g)),zorder=3)
        a.axvline(0,c="#AAA",lw=.5); a.set_xlim(-2,105); a.set_xlabel("Mutant frequency (%)"); a.grid(axis="x",color="#E7E7E7",lw=.5); panel(a,letter,facet); a.legend(frameon=False,loc="lower right",handletextpad=.3)
    axs[0].set_yticks(range(len(keys))); axs[0].set_yticklabels([MARKERS[k] for k in keys[::-1]])
    fig.suptitle("Figure 3 | Resistance profiles across ecology and vector species",x=.08,ha="left",fontsize=13,fontweight="bold")
    fig.text(.08,.015,"Points show mutant frequency; horizontal lines are Wilson 95% CIs; point size reflects the callable denominator.",fontsize=7,color="#444")
    fig.subplots_adjust(left=.15,right=.98,bottom=.14,top=.81,wspace=.12); export(fig,out,"figure3_ecological_abacus")

def figure4(d,coi,out):
    q=d.merge(coi[["sample_id","infection_class"]],on="sample_id",how="left")
    q["complexity"]=q.infection_class.map({"mono_compatible":"No heterozygosity detected","polygenomic":"Polygenomic signal"})
    q["callable_loci"]=q[[c for c in q if c.endswith("_coverage_above_threshold")]].fillna(False).sum(axis=1)
    src=q[["sample_id","sibling_species","host_feeding_type","pf_ct","infection_class","complexity","callable_loci"]]; src.to_csv(out/"figure4_complexity_source.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.2,4.4),gridspec_kw={"width_ratios":[1.25,1],"wspace":.35})
    ax=axs[0]; combos=[(s,h) for h in ["single_host","mixed_host"] for s in SPECIES_ORDER]; x=np.arange(len(combos)); bottom=np.zeros(len(combos));
    for cat,col in [("No heterozygosity detected",OI["sky"]),("Polygenomic signal",OI["vermillion"])]:
        vals=[]
        for s,h in combos:
            z=q[(q.sibling_species==s)&(q.host_feeding_type==h)].complexity.dropna(); vals.append((z==cat).mean() if len(z) else 0)
        ax.bar(x,100*np.array(vals),bottom=100*bottom,color=col,label=cat,width=.76); bottom+=vals
    short={"An_arabiensis":"Arabiensis","An_coluzzii":"Coluzzii","An_gambiae_s.s":"Gambiae s.s."}
    ax.set_xticks(x); ax.set_xticklabels([short[s]+"\n"+("single" if h=="single_host" else "mixed") for s,h in combos],rotation=0,ha="center",fontsize=6.7)
    ax.set_ylabel("Composition (%)"); ax.set_ylim(0,100); ax.legend(frameon=False,loc="upper left"); panel(ax,"A","Polygenomic signal by vector and feeding type")
    ax=axs[1]
    cats=["single_host","mixed_host"]; data=[q.loc[q.host_feeding_type==h,"pf_ct"].dropna().to_numpy() for h in cats]
    bp=ax.boxplot(data,positions=[0,1],widths=.5,patch_artist=True,showfliers=False,medianprops={"color":"black"})
    for b,c in zip(bp["boxes"],[OI["green"],OI["orange"]]): b.set_facecolor(c); b.set_alpha(.75)
    rng=np.random.default_rng(19)
    for i,z in enumerate(data): ax.scatter(i+rng.uniform(-.15,.15,len(z)),z,s=5,c="#333",alpha=.20,rasterized=True)
    ax.set_xticks([0,1]); ax.set_xticklabels([f"Single host\n$n$={len(data[0])}",f"Mixed host\n$n$={len(data[1])}"]); ax.set_ylabel(r"$P. falciparum$ qPCR Ct"); panel(ax,"B","Mixed-host dilution hypothesis")
    ax.text(.02,.02,"Higher Ct indicates lower parasite template.",transform=ax.transAxes,fontsize=7,color="#444")
    fig.suptitle("Figure 4 | Vector ecology, parasite complexity and template dilution",x=.06,ha="left",fontsize=13,fontweight="bold")
    fig.text(.06,.01,"Polygenomic signal = ≥1 non-artefactual heterozygous call; absence is not proof of monoclonality.",fontsize=7,color="#444")
    fig.subplots_adjust(left=.08,right=.98,bottom=.18,top=.78,wspace=.35); export(fig,out,"figure4_complexity_and_dilution")

def figure5(d,out):
    # A437G has enough mutant and wild-type calls for adjusted inference; K76T has zero events and is not modelled.
    d=d.copy(); d["A437G"]=(pd.to_numeric(d["dhps_1306_TG"],errors="coerce")>0).where(d["dhps_1306_TG"].notna())
    callcols=[c for c in d if c.endswith("_coverage_above_threshold")]; d["all_six_callable"]=(d[callcols].fillna(False).sum(axis=1)==len(callcols)).astype(int)
    a,u1=logistic_ridge(d,"A437G"); a["outcome"]=r"$pfdhps$ A437G carriage"
    b,u2=logistic_ridge(d,"all_six_callable"); b["outcome"]="All six loci callable"
    s=pd.concat([a,b],ignore_index=True); s.to_csv(out/"figure5_site_adjusted_models.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.2,4.6),sharey=False,gridspec_kw={"wspace":.10})
    for ax,(outcome,z),letter,col in zip(axs,s.groupby("outcome",sort=False),"AB",[OI["vermillion"],OI["blue"]]):
        z=z.reset_index(drop=True); y=np.arange(len(z))[::-1]
        lo=np.clip(z.conf_low,0.03,30); hi=np.clip(z.conf_high,0.03,30); ax.hlines(y,lo,hi,color=col,lw=1.4); ax.scatter(z.odds_ratio,y,s=30,c=col,edgecolor="white",lw=.5,zorder=3)
        ax.axvline(1,c="#555",ls="--",lw=.8); ax.set_xscale("log"); ax.set_xlim(.03,30); ax.set_xlabel("Adjusted odds ratio (log scale)"); ax.grid(axis="x",which="both",color="#E6E6E6",lw=.5); panel(ax,letter,outcome+f" ($n$={int(z.n.iloc[0])})")
        ax.set_yticks(y)
        if ax is axs[0]: ax.set_yticklabels(z.label)
        else: ax.set_yticklabels([]); ax.tick_params(axis="y",length=0)
    fig.suptitle("Figure 5 | Ecological predictors of resistance detection and assay success",x=.06,ha="left",fontsize=13,fontweight="bold")
    fig.text(.06,.015,"Site-adjusted penalised logistic models; points are ORs and lines approximate Wald 95% CIs. Reference: Forest, An. coluzzii, single-host feeding.",fontsize=7,color="#444")
    fig.subplots_adjust(left=.25,right=.98,bottom=.16,top=.78,wspace=.10); export(fig,out,"figure5_adjusted_effects")

def main():
    p=argparse.ArgumentParser(); p.add_argument("--result-dir",required=True,type=Path); p.add_argument("--metadata",required=True,type=Path); p.add_argument("--geo-dir",type=Path); p.add_argument("--outdir",required=True,type=Path); a=p.parse_args()
    style(); a.outdir.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(a.metadata); req={"sample_id","bioclimatic_zone","sibling_species","collection_site","latitude","longitude","host_feeding_type","pf_ct"}; missing=req-set(meta)
    if missing: raise SystemExit(f"Clean metadata missing required columns: {sorted(missing)}")
    # The clean file retains repeated run controls (KH2/NC) for laboratory
    # traceability. They have no ecological metadata and are not study units.
    if "sample_status" in meta:
        meta=meta[~meta.sample_status.isin(["PC","NC"])].copy()
    if meta.sample_id.duplicated().any(): raise SystemExit("Clean metadata sample_id must be unique")
    summary=a.result_dir/"04_summary"; variants=pd.read_csv(summary/"tables/summary_drug_resistance_variants_.csv"); cov=pd.read_csv(summary/"coverage/summary_coverage_by_run_sample.csv")
    d=variants.drop(columns=["collection_site"],errors="ignore").merge(meta,on="sample_id",how="left",validate="one_to_one").merge(cov.drop(columns=["run_name"],errors="ignore"),on="sample_id",how="left",validate="one_to_one")
    if d.bioclimatic_zone.isna().all(): raise SystemExit("No analyzed sample IDs matched clean metadata")
    coi_file=a.result_dir/"05_complexity/complexity_of_infection.tsv"; coi=pd.read_csv(coi_file,sep="\t")
    mapfile=summary/"processed_genotypes/summary_variants_linked_to_metadata_intermediate_file.csv"; mp=pd.read_csv(mapfile)[["sample_id","ont_multiplex_group","ont_barcode"]].drop_duplicates()
    coi=coi.merge(mp,left_on=["run_name","barcode"],right_on=["ont_multiplex_group","ont_barcode"],how="left").drop_duplicates("sample_id")
    adm0=load_geo(a.geo_dir/"ghana_ADM0.geojson" if a.geo_dir else None); adm1=load_geo(a.geo_dir/"ghana_ADM1.geojson" if a.geo_dir else None)
    figure1(d,adm0,adm1,a.outdir); figure2(d,a.outdir); figure3(d,a.outdir); figure4(d,coi,a.outdir); figure5(d,a.outdir)
    (a.outdir/"README.md").write_text("""# DRAG1 roadmap-driven publication figures\n\nFigures follow `Bioinformatics Visualization Roadmap Development.md` and use `samplesheet_clean_metadata.csv` exclusively for ecological metadata.\n\n1. Geographic pfdhfr haplotype composition.\n2. Clinically relevant pfdhfr/pfdhps mutation intersections.\n3. Resistance-frequency abacus by zone and vector species.\n4. Polygenomic signal and mixed-host dilution.\n5. Site-adjusted ecological effect sizes.\n\nAll mutation denominators are callable specimens only. Missing calls are never coded as wild type. PNG (400 dpi), SVG and PDF are supplied with source CSVs.\n""")
    print(f"Wrote five roadmap-driven publication figures to {a.outdir}")
if __name__=="__main__": main()

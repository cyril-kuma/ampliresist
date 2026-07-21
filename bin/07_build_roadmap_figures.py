#!/usr/bin/env python3
"""Roadmap-driven publication figures for DRAG1 molecular xenomonitoring.

The figures deliberately follow the narrative in
``Bioinformatics Visualization Roadmap Development.md``.  Metadata are read
only from the cleaned CSV supplied with ``--metadata``.  Genomic calls and
callability remain products of the pipeline and are never inferred from
missing VCF records.

DENOMINATOR CONTRACT
--------------------
Every panel that reports a proportion also reports the denominator it used AND
the cohort total it was drawn from, so a reader can always see how much of the
cohort a number rests on.  Three rules make that hold:

  * ``n_total`` is never recomputed here.  It is read from stage 3's frequency
    table, which is the pipeline's authority on the analysed cohort size.
  * A marker with zero mutant calls is a RESULT (0% of n callable), not a
    missing marker.  Markers are dropped from a figure only when they are
    ``not_callable`` — and those are then named explicitly rather than vanishing.
  * Marker labels and reference polarity are read from stage 2's
    ``*_marker_polarity.csv``.  The wide variants table counts NON-REFERENCE
    alleles; for markers where 3D7 itself carries the mutant (dhps A437G) that
    is the inverse of carriage, so it is flipped here exactly as stage 3 does.
"""
from __future__ import annotations

import argparse, json, math, re, textwrap
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Rectangle
from matplotlib.lines import Line2D

OI = {"orange":"#E69F00","sky":"#56B4E9","green":"#009E73","yellow":"#F0E442",
      "blue":"#0072B2","vermillion":"#D55E00","purple":"#CC79A7","black":"#111111"}
GREY = "#BDBDBD"
ZONE_ORDER = ["Coastal_Savannah","Forest","Northern_Savannah"]
ZONE_LABEL = {"Coastal_Savannah":"Coastal Savannah","Forest":"Forest","Northern_Savannah":"Northern Savannah"}
SPECIES_ORDER = ["An_arabiensis","An_coluzzii","An_gambiae_s.s"]
SPECIES_LABEL = {"An_arabiensis":r"$An.$ arabiensis","An_coluzzii":r"$An.$ coluzzii","An_gambiae_s.s":r"$An.$ gambiae s.s."}
GENES = ["crt","dhfr","dhps","mdr1","k13","csp"]
# Display order only. Membership is decided by the pipeline's marker tables, so a
# marker added upstream still appears here without editing this file.
MARKER_ORDER = ["crt_227_AC","dhfr_152_AT","dhfr_175_TC","dhfr_323_GA","dhfr_323_GC",
                "dhfr_490_AT","dhfr_492_AG","dhps_1306_TG","dhps_1307_CT","dhps_1310_GC",
                "dhps_1618_AG","dhps_1620_AT","dhps_1742_CG","dhps_1837_GA","dhps_1837_GT",
                "mdr1_256_AT","mdr1_551_AT"]

def style():
    mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"axes.titlesize":10,
        "axes.labelsize":9,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5,
        "axes.linewidth":0.7,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
        "savefig.dpi":400,"savefig.bbox":"tight","axes.spines.top":False,"axes.spines.right":False})

def panel(ax, letter, title, size=None):
    ax.text(-.10,1.015,letter,transform=ax.transAxes,fontsize=12,fontweight="bold",va="bottom")
    ax.set_title(title,loc="left",fontweight="bold",pad=5,fontsize=size)

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

# ---------------------------------------------------------------------------
# Marker identity, polarity and cohort denominators — all read from the pipeline
# ---------------------------------------------------------------------------

def pretty_marker(aa_mut_name):
    """'DHPS_A437G' -> '$pfdhps$ A437G'."""
    gene,_,mut = str(aa_mut_name).partition("_")
    return rf"$pf{gene.lower()}$ {mut}" if mut else str(aa_mut_name)

def load_polarity(genotype_dir):
    """Stage 2's marker table: canonical amino-acid names + reference polarity.

    Hard-coding these labels is how the old figure suite came to print three dhps
    markers one row out of step with the catalogue (S436A shown as A437G, A437G
    shown as K540E, A613S shown as A581G).  They are derived now.
    """
    hits=sorted(Path(genotype_dir).glob("*_marker_polarity.csv"))
    if not hits: raise SystemExit(f"No *_marker_polarity.csv in {genotype_dir}; cannot label markers or correct polarity")
    p=pd.read_csv(hits[0])
    for col in ("snp_id","aa_mut_name","ref_is_mutant"):
        if col not in p: raise SystemExit(f"{hits[0].name} is missing column '{col}'")
    p["label"]=p.aa_mut_name.map(pretty_marker)
    p["ref_is_mutant"]=p.ref_is_mutant.astype(str).str.upper().isin(["TRUE","T","1"])
    return p

def carriage(d, polarity):
    """Polarity-corrected carriage of each NAMED mutation.

    The wide variants table stores counts of NON-REFERENCE alleles.  Where 3D7
    carries the mutant allele (dhps A437G: reference codon GGT = Gly = 437G), a
    non-reference call is the WILD TYPE, so 'raw > 0' is carriage inverted.
    Stage 3 corrects this as n_callable - Nref_count; the same flip is applied
    here so the figures cannot disagree with the frequency table.
    1 = carries mutation, 0 = does not, NaN = not callable.
    """
    flip=set(polarity.loc[polarity.ref_is_mutant,"snp_id"])
    out={}
    for snp in polarity.snp_id:
        if snp not in d: continue
        raw=pd.to_numeric(d[snp],errors="coerce")
        val=(raw>0).astype(float).where(raw.notna())
        out[snp]=(1.0-val) if snp in flip else val
    return pd.DataFrame(out,index=d.index)

def load_frequencies(summary_dir):
    """Stage 3's frequency table — the authority on n_callable and n_total."""
    f=summary_dir/"tables/summary_drug_resistance_frequencies_.csv"
    if not f.exists(): raise SystemExit(f"Missing frequency table {f}; it carries the cohort denominator n_total")
    fr=pd.read_csv(f)
    for col in ("SNP","n_callable","n_total","callable"):
        if col not in fr: raise SystemExit(f"Frequency table is missing column '{col}'")
    return fr

def order_markers(snps, polarity):
    known=[s for s in MARKER_ORDER if s in snps]
    return known+[s for s in polarity.snp_id if s in snps and s not in known]

def wrap(s, width=128):
    """Hard-wrap footnotes.

    savefig uses bbox='tight', so a single long fig.text line silently widens the
    exported canvas to fit it — a 7.2 in figure came out 19 in wide and squashed.
    """
    return "\n".join(textwrap.fill(ln,width) for ln in s.split("\n"))

def denom_note(n_used, n_total, what="specimens"):
    pct=100*n_used/n_total if n_total else float("nan")
    return f"$n$={n_used:,} of {n_total:,} {what} ({pct:.1f}%)"

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def figure0(d, car, freq, n_total, min_cov, out):
    """Cohort accounting: what every other denominator in the suite is drawn from.

    Without this panel the suite reported denominators between 106 and 438 across
    five figures with no way to tell which numbers shared one, or how much of the
    cohort any of them covered.
    """
    callcols=[f"{g}_coverage_above_threshold" for g in GENES if f"{g}_coverage_above_threshold" in d]
    qc_pass=int((d.sample_qc=="PASS").sum()) if "sample_qc" in d else n_total
    sp=[k for k in ["dhfr_152_AT","dhfr_175_TC","dhfr_323_GA","dhps_1306_TG","dhps_1310_GC"] if k in car]

    # The funnel is chosen to explain the denominators the other figures use:
    # figure 2 is the complete-case row, figure 5's assay-success model the last.
    flow=[("Specimens in\nresistance cohort",n_total),
          ("Sequence data\nrecovered",int(d[callcols].notna().any(axis=1).sum()) if callcols else n_total),
          ("Sample QC\nPASS",qc_pass),
          ("SP complete case\n(Fig. 2)",int(car[sp].notna().all(axis=1).sum()) if sp else 0),
          ("All 6 amplicons\ncallable (Fig. 5)",int((d[callcols].fillna(False).astype(bool).sum(axis=1)==len(callcols)).sum()) if callcols else 0)]

    gene_rows=[]
    for g in GENES:
        c=f"{g}_coverage_above_threshold"
        if c in d: gene_rows.append({"gene":g,"callable":int(d[c].fillna(False).astype(bool).sum()),"n_total":n_total})
    gene_df=pd.DataFrame(gene_rows)

    src=pd.DataFrame([{"stage":s,"n":n,"n_total":n_total,"pct":100*n/n_total} for s,n in flow])
    pd.concat([src.assign(block="cohort_flow"),
               gene_df.assign(block="amplicon_callability",stage=gene_df.gene,n=gene_df.callable,
                              pct=100*gene_df.callable/n_total)[["block","stage","n","n_total","pct"]],
               freq.assign(block="marker_denominator",stage=freq.SNP,n=freq.Nref_count,
                           pct=100*freq.Nref_count/freq.n_callable)[["block","stage","n","n_total","pct"]]],
              ignore_index=True).to_csv(out/"figure0_cohort_accounting_source.csv",index=False)

    fig=plt.figure(figsize=(7.2,8.2)); gs=fig.add_gridspec(3,1,height_ratios=[1.0,0.85,1.85],hspace=.42)

    ax=fig.add_subplot(gs[0]); x=np.arange(len(flow)); vals=[n for _,n in flow]
    ax.bar(x,vals,color=[OI["blue"]]+[OI["sky"]]*2+[OI["green"],OI["orange"]],width=.68)
    ax.axhline(n_total,color="#555",ls="--",lw=.8)
    ax.text(len(flow)-.45,n_total*1.03,f"cohort $n$={n_total:,}",va="bottom",ha="right",fontsize=7,color="#555")
    for i,v in enumerate(vals):
        ax.text(i,v+n_total*.02,f"{v:,}\n{100*v/n_total:.1f}%",ha="center",va="bottom",fontsize=7,fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([s for s,_ in flow],fontsize=6.6); ax.set_ylabel("Specimens")
    ax.set_ylim(0,n_total*1.28); panel(ax,"A","Cohort accounting — every specimen analysed")

    ax=fig.add_subplot(gs[1]); x=np.arange(len(gene_df))
    ax.bar(x,gene_df.callable,color=OI["sky"],width=.68)
    ax.bar(x,gene_df.n_total-gene_df.callable,bottom=gene_df.callable,color=GREY,width=.68)
    for i,r in gene_df.iterrows():
        ax.text(i,r.n_total*1.03,f"{int(r.callable)}/{int(r.n_total)}",ha="center",va="bottom",fontsize=7,fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([f"$pf{g}$" for g in gene_df.gene],fontsize=8)
    ax.set_ylabel("Specimens"); ax.set_ylim(0,n_total*1.22)
    panel(ax,"B",f"Amplicon callability at ≥{min_cov:g}× median — why denominators differ")

    ax=fig.add_subplot(gs[2])
    fr=freq.set_index("SNP").reindex(order_markers(list(freq.SNP),POLARITY)).reset_index()
    y=np.arange(len(fr))[::-1]
    ax.barh(y,fr.n_callable.fillna(0),color=OI["sky"],height=.66)
    ax.barh(y,fr.n_total-fr.n_callable.fillna(0),left=fr.n_callable.fillna(0),color=GREY,height=.66)
    ax.barh(y,fr.Nref_count.fillna(0),color=OI["vermillion"],height=.30)
    for yi,(_,r) in zip(y,fr.iterrows()):
        ok=r.callable=="callable"
        ax.text(n_total*1.03,yi,f"{int(r.Nref_count or 0)}/{int(r.n_callable or 0)}" if ok else "not callable",
                va="center",fontsize=6.6,fontweight="bold",color="#333" if ok else OI["vermillion"])
    ax.set_yticks(y); ax.set_yticklabels([MARKER_LABEL.get(s,s) for s in fr.SNP],fontsize=7)
    ax.set_xlim(0,n_total*1.30); ax.set_xlabel("Specimens")
    panel(ax,"C","Per-marker denominators (mutant / callable, of cohort)")

    # One legend for panels B and C, which share an encoding, placed clear of both.
    fig.legend([Rectangle((0,0),1,1,fc=OI["sky"]),Rectangle((0,0),1,1,fc=GREY),Rectangle((0,0),1,1,fc=OI["vermillion"])],
               ["Callable","Not callable","Carrying mutation"],
               frameon=False,ncol=3,fontsize=7.5,loc="lower center",bbox_to_anchor=(.55,.012))
    fig.suptitle("Figure 0 | How much of the cohort each result rests on",x=.06,ha="left",fontsize=13,fontweight="bold")
    fig.text(.06,-.01,wrap("Grey = not callable, never counted as wild type. A marker with 0 mutants among callable specimens is a result (0%), not a missing marker.\n"
                      r"$pfdhps$ A437G is reference-polarity corrected: 3D7 itself carries 437G, so a non-reference call there is the wild type."),fontsize=6.8,color="#444",va="top")
    fig.subplots_adjust(left=.21,right=.95,bottom=.085,top=.93); export(fig,out,"figure0_cohort_accounting")

def figure1(d, car, n_total, adm0, adm1, out):
    # Site-level DHFR haplotypes communicate geography and clinically meaningful composition.
    # Haplotypes are taken from the data: the old hard-coded list of four silently
    # discarded every NCNI specimen (20 of 126 callable) from the map.
    obs=d.dhfr_haplotype.dropna()
    haps=[h for h in ["NCSI","NCNI","NRNI","IRNI","ICNI","IRSI","ICSI"] if h in set(obs)]
    haps+=sorted(set(obs)-set(haps))
    pal=[GREY,OI["yellow"],OI["sky"],OI["orange"],OI["purple"],OI["green"],OI["vermillion"],OI["blue"]]
    hc={h:pal[i%len(pal)] for i,h in enumerate(haps)}
    q=d.dropna(subset=["dhfr_haplotype","collection_site","latitude","longitude"])
    site_tot=d.dropna(subset=["collection_site"]).groupby("collection_site").size().rename("n_site")
    ct=q.groupby(["collection_site","latitude","longitude","bioclimatic_zone","dhfr_haplotype"]).size().unstack(fill_value=0).reindex(columns=haps,fill_value=0).reset_index()
    ct["n_callable"]=ct[haps].sum(axis=1)
    ct=ct.merge(site_tot,on="collection_site",how="left")
    ct["n_total_cohort"]=n_total; ct.to_csv(out/"figure1_geographic_source.csv",index=False)
    fig,ax=plt.subplots(figsize=(7.2,6.2)); panel(ax,"A","Geographic composition of pfdhfr haplotypes")
    for ring in adm0: ax.fill([p[0] for p in ring],[p[1] for p in ring],fc="#FAFAFA",ec="#555",lw=.8,zorder=0)
    for ring in adm1: ax.plot([p[0] for p in ring],[p[1] for p in ring],c="#D0D0D0",lw=.45,zorder=.5)
    nudges={"Aplaku":(-.18,-.05),"Teshie":(.18,.04),"Obuasi":(-.12,.04),"Bekwai":(.12,.08),
            "Chiraa":(-.30,.14),"Sunyani":(.24,-.14),"Berekum":(-.34,-.06)}
    for _,r in ct.sort_values("n_callable").iterrows():
        x,y=float(r.longitude),float(r.latitude); dx,dy=nudges.get(r.collection_site,(0,0)); x2,y2=x+dx,y+dy
        n=int(r.n_callable); tot=int(r.n_site); rad=.035+.018*math.sqrt(n)
        if dx or dy: ax.plot([x,x2],[y,y2],c="#777",lw=.5,zorder=2)
        start=90
        for h in haps:
            extent=360*int(r[h])/n
            ax.add_patch(Wedge((x2,y2),rad,start,start+extent,facecolor=hc[h],edgecolor="white",lw=.5,zorder=3)); start+=extent
        ax.add_patch(Circle((x2,y2),rad,fill=False,ec="#333",lw=.5,zorder=4))
        ax.annotate(f"{r.collection_site}\n$n$={n}/{tot}",(x2,y2),xytext=(rad*1.15,0),textcoords="offset points",fontsize=6.5,va="center",fontweight="bold")
    # Sites with specimens but no dhfr-callable specimen would otherwise be absent
    # from the map entirely, which reads as "not sampled" rather than "not callable".
    silent=sorted(set(site_tot.index)-set(ct.collection_site))
    if silent:
        sl=d.dropna(subset=["collection_site","latitude","longitude"]).drop_duplicates("collection_site").set_index("collection_site")
        for s in silent:
            if s not in sl.index: continue
            x,y=float(sl.loc[s,"longitude"]),float(sl.loc[s,"latitude"]); dx,dy=nudges.get(s,(0,0))
            ax.scatter(x,y,marker="x",s=26,c="#B00020",lw=1.1,zorder=5)
            if dx or dy: ax.plot([x,x+dx],[y,y+dy],c="#B00020",lw=.4,zorder=2)
            ax.annotate(f"{s}\n0/{int(site_tot[s])}",(x+dx,y+dy),xytext=(4,0),textcoords="offset points",
                        fontsize=6,va="center",color="#B00020")
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)"); ax.set_aspect("equal"); ax.grid(False)
    handles=[Rectangle((0,0),1,1,fc=hc[h]) for h in haps]+[Line2D([],[],marker="x",ls="none",c="#B00020")]
    ax.legend(handles,haps+["sampled, none callable"],title=r"$pfdhfr$ haplotype",frameon=False,loc="upper left",bbox_to_anchor=(1.01,.98))
    ax.text(.01,.01,f"Pie area ∝ DHFR-callable specimens; slices show within-site composition.\nLabels are callable/sampled per site. Cohort: {denom_note(len(q),n_total)} are $pfdhfr$-callable.\nSmall denominators are descriptive, not population prevalence.",transform=ax.transAxes,fontsize=7,color="#444")
    fig.suptitle("Figure 1 | Where are pyrimethamine-resistance haplotypes circulating?",x=.08,ha="left",fontsize=13,fontweight="bold")
    fig.subplots_adjust(left=.12,right=.78,bottom=.11,top=.86); export(fig,out,"figure1_geographic_resistance")

def figure2(d, car, n_total, out):
    # Intersections are inherently complete-case; the denominator is stated on the
    # figure so the reader can see how much of the cohort the pattern rests on.
    keys=[k for k in ["dhfr_152_AT","dhfr_175_TC","dhfr_323_GA","dhps_1306_TG","dhps_1310_GC"] if k in car]
    q=car.dropna(subset=keys).copy(); counts=q[keys].astype(int).astype(str).agg("".join,axis=1).value_counts().head(15)
    src=[]
    for pattern,n in counts.items():
        row={"pattern":pattern,"n":n,"n_complete_case":len(q),"n_total":n_total}
        row.update({k:int(v) for k,v in zip(keys,pattern)}); src.append(row)
    src=pd.DataFrame(src); src.to_csv(out/"figure2_upset_source.csv",index=False)
    fig=plt.figure(figsize=(7.2,5.6)); gs=fig.add_gridspec(2,1,height_ratios=[2.2,1],hspace=.04); ax=fig.add_subplot(gs[0]); mx=fig.add_subplot(gs[1],sharex=ax)
    x=np.arange(len(src)); ax.bar(x,src.n,color=OI["blue"],width=.75); ax.set_ylabel("Specimens"); ax.set_xticks([]); panel(ax,"A","Most frequent clinically relevant mutation intersections")
    for i,n in enumerate(src.n): ax.text(i,n+max(src.n)*.015,f"{n}\n{100*n/len(q):.0f}%",ha="center",va="bottom",fontsize=6.6,fontweight="bold")
    for row,k in enumerate(keys):
        yy=len(keys)-1-row; mx.text(-.75,yy,MARKER_LABEL.get(k,k),ha="right",va="center",fontsize=8)
        mx.scatter(x,np.full(len(x),yy),s=15,c="#D8D8D8",zorder=1)
        on=np.where(src[k].to_numpy()==1)[0]; mx.scatter(on,np.full(len(on),yy),s=24,c=OI["black"],zorder=2)
        for i in range(len(src)):
            ys=np.where(src.loc[i,keys].to_numpy()==1)[0]
            if len(ys)>1: mx.plot([i,i],[len(keys)-1-ys.max(),len(keys)-1-ys.min()],c=OI["black"],lw=1.2,zorder=1)
    mx.set_ylim(-.7,len(keys)-.3); mx.set_yticks([]); mx.set_xticks(x); mx.set_xticklabels([f"I{i+1}" for i in x],fontsize=7); mx.set_xlabel("Mutation intersection (complete calls only)"); mx.spines[["left","right","top"]].set_visible(False)
    fig.suptitle("Figure 2 | Multigenic architecture of SP resistance",x=.08,ha="left",fontsize=13,fontweight="bold")
    lim=min((int(car[k].notna().sum()) for k in keys),default=0)
    fig.text(.08,-.01,wrap(f"Complete-case denominator: {denom_note(len(q),n_total)}; percentages are of the complete-case set. The binding constraint is $pfdhfr$ callability ({lim}/{n_total}). Grey = absent; black = present. Missing calls are excluded, never treated as wild type."),fontsize=6.8,color="#444",va="top")
    fig.tight_layout(rect=(.12,.02,1,.94)); export(fig,out,"figure2_multigenic_upset")

def figure3(d, car, freq, n_total, out):
    # Every CALLABLE marker is plotted. The previous filter kept only markers whose
    # mutant count was > 0, which deleted pfcrt K76T, pfmdr1 N86Y and four dhps
    # markers — all with 413-447 callable specimens — because their result was 0%.
    fr=freq.set_index("SNP")
    keys=[k for k in order_markers([c for c in car if car[c].notna().any()],POLARITY)]
    not_callable=[MARKER_LABEL.get(s,s) for s in fr.index[fr.callable.ne("callable")]]
    rows=[]
    for facet,var,levels in [("Bioclimatic zone","bioclimatic_zone",ZONE_ORDER),("Vector species","sibling_species",SPECIES_ORDER)]:
        for level in levels:
            grp=d[var]==level; n_grp=int(grp.sum())
            for k in keys:
                z=car.loc[grp,k].dropna(); n=len(z); mut=int(z.sum()); lo,hi=wilson(mut,n)
                rows.append({"facet":facet,"group":level,"marker":k,"label":MARKER_LABEL.get(k,k),
                             "mutant":mut,"callable":n,"n_group":n_grp,"n_total":n_total,
                             "frequency":mut/n if n else np.nan,"low":lo,"high":hi})
    s=pd.DataFrame(rows); s.to_csv(out/"figure3_abacus_source.csv",index=False)
    h=max(4.8,1.6+0.30*len(keys))
    fig,axs=plt.subplots(1,2,figsize=(7.2,h),sharey=True,gridspec_kw={"wspace":.12})
    for a,(facet,sub),letter in zip(axs,s.groupby("facet",sort=False),"AB"):
        groups=list(dict.fromkeys(sub.group)); ypos={k:i for i,k in enumerate(keys[::-1])}
        offsets=np.linspace(-.26,.26,len(groups)); colors=[OI["sky"],OI["green"],OI["orange"]]
        for off,g,c in zip(offsets,groups,colors):
            z=sub[sub.group==g]
            yy=np.array([ypos[k] for k in z.marker])+off; xx=100*z.frequency.to_numpy()
            a.hlines(yy,100*z.low,100*z.high,color=c,lw=1)
            a.scatter(xx,yy,s=14+2*np.sqrt(z.callable),c=c,edgecolor="white",lw=.35,
                      label=f"{ZONE_LABEL.get(g,SPECIES_LABEL.get(g,g))} ($n$={int(z.n_group.iloc[0])})",zorder=3)
        a.axvline(0,c="#AAA",lw=.5); a.set_xlim(-3,106); a.set_xlabel("Mutant frequency (%)")
        a.grid(axis="x",color="#E7E7E7",lw=.5); panel(a,letter,facet)
        a.legend(frameon=False,loc="lower right",handletextpad=.3,fontsize=6.6)
    # Cohort-wide callable denominator against each marker: the reader can see that
    # a 0% at n=413 and a 44% at n=126 are not equally well supported.
    axs[0].set_yticks(range(len(keys)))
    axs[0].set_yticklabels([f"{MARKER_LABEL.get(k,k)}\n$n$={int(fr.n_callable.get(k,0) or 0)}/{n_total}" for k in keys[::-1]],fontsize=6.8)
    fig.suptitle("Figure 3 | Resistance profiles across ecology and vector species",x=.08,ha="left",fontsize=13,fontweight="bold")
    foot=("Points show mutant frequency; horizontal lines are Wilson 95% CIs; point size reflects the callable denominator. "
          "Markers at 0% are results, not missing data. $pfdhps$ A437G is reference-polarity corrected (3D7 carries 437G).")
    if not_callable: foot+=" Not callable anywhere, so not plotted: "+", ".join(not_callable)+"."
    fig.text(.08,-.01,wrap(foot),fontsize=6.6,color="#444",va="top")
    fig.subplots_adjust(left=.22,right=.98,bottom=.09,top=.85,wspace=.12); export(fig,out,"figure3_ecological_abacus")

def figure4(d,coi,n_total,out):
    q=d.merge(coi[["sample_id","infection_class"]],on="sample_id",how="left")
    q["complexity"]=q.infection_class.map({"mono_compatible":"No heterozygosity detected","polygenomic":"Polygenomic signal"})
    q["callable_loci"]=q[[c for c in q if c.endswith("_coverage_above_threshold")]].fillna(False).astype(bool).sum(axis=1)
    src=q[["sample_id","sibling_species","host_feeding_type","pf_ct","infection_class","complexity","callable_loci"]].copy()
    src["n_total"]=n_total; src.to_csv(out/"figure4_complexity_source.csv",index=False)
    fig,axs=plt.subplots(1,2,figsize=(7.2,4.6),gridspec_kw={"width_ratios":[1.25,1],"wspace":.35})
    ax=axs[0]; combos=[(s,h) for h in ["single_host","mixed_host"] for s in SPECIES_ORDER]; x=np.arange(len(combos)); bottom=np.zeros(len(combos))
    ns=[int(q[(q.sibling_species==s)&(q.host_feeding_type==h)].complexity.notna().sum()) for s,h in combos]
    for cat,col in [("No heterozygosity detected",OI["sky"]),("Polygenomic signal",OI["vermillion"])]:
        vals=[]
        for s,h in combos:
            z=q[(q.sibling_species==s)&(q.host_feeding_type==h)].complexity.dropna(); vals.append((z==cat).mean() if len(z) else 0)
        ax.bar(x,100*np.array(vals),bottom=100*bottom,color=col,label=cat,width=.76); bottom+=vals
    short={"An_arabiensis":"Arabiensis","An_coluzzii":"Coluzzii","An_gambiae_s.s":"Gambiae s.s."}
    # Stacked percentages with no denominator invite reading a 2-specimen cell as a
    # 100-specimen one; n is printed on every bar.
    for i,n in enumerate(ns): ax.text(i,101,f"$n$={n}",ha="center",va="bottom",fontsize=6.4,fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([short[s]+"\n"+("single" if h=="single_host" else "mixed") for s,h in combos],rotation=0,ha="center",fontsize=6.7)
    ax.set_ylabel("Composition (%)"); ax.set_ylim(0,112); ax.set_yticks([0,25,50,75,100])
    ax.legend(frameon=False,loc="lower left",bbox_to_anchor=(0,-.34),ncol=2,fontsize=7); panel(ax,"A","Polygenomic signal by vector and feeding type")
    ax=axs[1]
    cats=["single_host","mixed_host"]; data=[q.loc[q.host_feeding_type==h,"pf_ct"].dropna().to_numpy() for h in cats]
    bp=ax.boxplot(data,positions=[0,1],widths=.5,patch_artist=True,showfliers=False,medianprops={"color":"black"})
    for b,c in zip(bp["boxes"],[OI["green"],OI["orange"]]): b.set_facecolor(c); b.set_alpha(.75)
    rng=np.random.default_rng(19)
    for i,z in enumerate(data): ax.scatter(i+rng.uniform(-.15,.15,len(z)),z,s=5,c="#333",alpha=.20,rasterized=True)
    ax.set_xticks([0,1]); ax.set_xticklabels([f"Single host\n$n$={len(data[0])}",f"Mixed host\n$n$={len(data[1])}"]); ax.set_ylabel(r"$P. falciparum$ qPCR Ct"); panel(ax,"B","Mixed-host dilution hypothesis")
    ax.text(.02,.02,"Higher Ct indicates lower parasite template.",transform=ax.transAxes,fontsize=7,color="#444")
    fig.suptitle("Figure 4 | Vector ecology, parasite complexity and template dilution",x=.06,ha="left",fontsize=13,fontweight="bold")
    n_cx=int(q.complexity.notna().sum()); n_ct=int(q.pf_ct.notna().sum())
    fig.text(.06,-.01,wrap(f"Complexity classified for {denom_note(n_cx,n_total)}; Ct available for {denom_note(n_ct,n_total)}. Polygenomic signal = ≥1 non-artefactual heterozygous call; absence is not proof of monoclonality."),fontsize=6.8,color="#444",va="top")
    fig.subplots_adjust(left=.08,right=.98,bottom=.24,top=.78,wspace=.35); export(fig,out,"figure4_complexity_and_dilution")

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
        # Complete/quasi-complete separation: the likelihood has no interior
        # maximum and IRLS walks the coefficient to the clip bound, yielding an
        # odds ratio of ~5e8 with a standard error in the hundreds. Clipping that
        # to the axis renders a non-result as a maximal effect, so it is flagged
        # and drawn as not estimable instead.
        sep=bool(se[i]>=10 or abs(beta[i])>=10)
        rows.append({"term":nm,"label":labels.get(nm,nm),"log_or":beta[i],"se":se[i],"separated":sep,
                     "odds_ratio":np.nan if sep else math.exp(np.clip(beta[i],-20,20)),
                     "conf_low":np.nan if sep else math.exp(np.clip(beta[i]-1.96*se[i],-20,20)),
                     "conf_high":np.nan if sep else math.exp(np.clip(beta[i]+1.96*se[i],-20,20)),
                     "n":len(use),"events":int(y.sum())})
    return pd.DataFrame(rows), use

def figure5(d, car, freq, n_total, out):
    # An outcome with no variation cannot be modelled. The previous version fitted
    # dhps_1306_TG with 0 events in 405 specimens and drew six odds ratios of
    # exactly 1.0 with standard errors in the hundreds — an empty result rendered
    # as a finished figure. Outcomes are screened for variation first.
    d=d.copy()
    callcols=[c for c in d if c.endswith("_coverage_above_threshold")]
    d["all_six_callable"]=(d[callcols].fillna(False).astype(bool).sum(axis=1)==len(callcols)).astype(int)

    MIN_EVENTS=10
    cand=[]
    for k in order_markers([c for c in car if car[c].notna().any()],POLARITY):
        v=car[k]; ev=int(v.sum()); nn=int(v.notna().sum())
        if ev>=MIN_EVENTS and (nn-ev)>=MIN_EVENTS: cand.append((min(ev,nn-ev),k,ev,nn))
    skipped=[(MARKER_LABEL.get(k,k),int(car[k].sum()),int(car[k].notna().sum()))
             for k in car if car[k].notna().any() and not any(c[1]==k for c in cand)]
    blocks=[]
    if cand:
        _,k,ev,nn=max(cand)
        d["_outcome"]=car[k]
        a,_=logistic_ridge(d,"_outcome"); a["outcome"]=f"{MARKER_LABEL.get(k,k)} carriage"; blocks.append(a)
    b,_=logistic_ridge(d,"all_six_callable"); b["outcome"]="All six loci callable"; blocks.append(b)
    s=pd.concat(blocks,ignore_index=True); s["n_total"]=n_total
    s.to_csv(out/"figure5_site_adjusted_models.csv",index=False)
    fig,axs=plt.subplots(1,len(blocks),figsize=(7.2,4.6),sharey=False,squeeze=False,gridspec_kw={"wspace":.10})
    axs=axs[0]
    for ax,(outcome,z),letter,col in zip(axs,s.groupby("outcome",sort=False),"AB",[OI["vermillion"],OI["blue"]]):
        z=z.reset_index(drop=True); y=np.arange(len(z))[::-1]
        ok=~z.separated.to_numpy(); yo=y[ok]
        lo=np.clip(z.conf_low[ok],0.03,30); hi=np.clip(z.conf_high[ok],0.03,30)
        ax.hlines(yo,lo,hi,color=col,lw=1.4); ax.scatter(z.odds_ratio[ok],yo,s=30,c=col,edgecolor="white",lw=.5,zorder=3)
        for yi in y[~ok]:
            ax.hlines(yi,.03,30,color="#CCC",lw=.8,ls=":")
            ax.text(1,yi,"not estimable (separation)",ha="center",va="center",fontsize=6,color="#B00020",
                    bbox=dict(fc="white",ec="none",pad=1.2))
        ax.axvline(1,c="#555",ls="--",lw=.8); ax.set_xscale("log"); ax.set_xlim(.03,30); ax.set_xlabel("Adjusted odds ratio (log scale)"); ax.grid(axis="x",which="both",color="#E6E6E6",lw=.5)
        panel(ax,letter,outcome+f"\n{int(z.events.iloc[0])} events / $n$={int(z.n.iloc[0])} modelled (of {n_total})",size=8.5)
        ax.set_yticks(y)
        if ax is axs[0]: ax.set_yticklabels(z.label)
        else: ax.set_yticklabels([]); ax.tick_params(axis="y",length=0)
    fig.suptitle("Figure 5 | Ecological predictors of resistance detection and assay success",x=.06,ha="left",fontsize=13,fontweight="bold")
    foot=("Site-adjusted penalised logistic models; points are ORs and lines approximate Wald 95% CIs. Reference: Forest, An. coluzzii, single-host feeding.")
    if bool(s.separated.any()):
        foot+=("\nTerms marked not estimable show complete separation (every specimen in that stratum shares the outcome); "
               "no finite odds ratio exists and none is drawn.")
    if skipped:
        inv=", ".join(f"{lab} ({ev}/{nn})" for lab,ev,nn in skipped[:8])
        foot+=f"\nNot modelled — fewer than {MIN_EVENTS} events or non-events, so no variation to explain: {inv}."
    fig.text(.06,-.01,wrap(foot),fontsize=6.6,color="#444",va="top")
    fig.subplots_adjust(left=.25,right=.98,bottom=.13,top=.74,wspace=.10); export(fig,out,"figure5_adjusted_effects")

def main():
    p=argparse.ArgumentParser(); p.add_argument("--result-dir",required=True,type=Path); p.add_argument("--metadata",required=True,type=Path); p.add_argument("--geo-dir",type=Path); p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--min-cov",type=float,default=50,help="amplicon median-coverage threshold used upstream; annotated on figure 0")
    a=p.parse_args()
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

    global POLARITY, MARKER_LABEL
    POLARITY=load_polarity(a.result_dir/"02_genotype_calls")
    MARKER_LABEL=dict(zip(POLARITY.snp_id,POLARITY.label))
    freq=load_frequencies(summary)
    car=carriage(d,POLARITY)
    # n_total is the pipeline's own cohort size (stage 3), not len(d): if the two
    # ever disagree the figures are being built from a partial table and must stop.
    n_total=int(freq.n_total.max())
    if len(d)!=n_total:
        raise SystemExit(f"Cohort mismatch: {len(d)} specimens in the variants table but stage 3 reports n_total={n_total}. Figures would understate the denominator.")

    coi_file=a.result_dir/"05_complexity/complexity_of_infection.tsv"; coi=pd.read_csv(coi_file,sep="\t")
    mapfile=summary/"processed_genotypes/summary_variants_linked_to_metadata_intermediate_file.csv"; mp=pd.read_csv(mapfile)[["sample_id","ont_multiplex_group","ont_barcode"]].drop_duplicates()
    coi=coi.merge(mp,left_on=["run_name","barcode"],right_on=["ont_multiplex_group","ont_barcode"],how="left").drop_duplicates("sample_id")
    adm0=load_geo(a.geo_dir/"ghana_ADM0.geojson" if a.geo_dir else None); adm1=load_geo(a.geo_dir/"ghana_ADM1.geojson" if a.geo_dir else None)
    figure0(d,car,freq,n_total,a.min_cov,a.outdir)
    figure1(d,car,n_total,adm0,adm1,a.outdir); figure2(d,car,n_total,a.outdir); figure3(d,car,freq,n_total,a.outdir)
    figure4(d,coi,n_total,a.outdir); figure5(d,car,freq,n_total,a.outdir)
    (a.outdir/"README.md").write_text(f"""# DRAG1 roadmap-driven publication figures

Figures follow `Bioinformatics Visualization Roadmap Development.md` and use
`samplesheet_clean_metadata.csv` exclusively for ecological metadata.

0. Cohort accounting and per-marker denominators.
1. Geographic pfdhfr haplotype composition.
2. Clinically relevant pfdhfr/pfdhps mutation intersections.
3. Resistance-frequency abacus by zone and vector species.
4. Polygenomic signal and mixed-host dilution.
5. Site-adjusted ecological effect sizes.

## Denominators

Cohort size is n = {n_total} specimens, taken from stage 3's frequency table and
verified against the variants table before any figure is drawn. Every panel that
reports a proportion also reports the denominator it used and the cohort total.

Mutation denominators are callable specimens only; missing calls are never coded
as wild type. Markers with zero mutant calls are reported at 0% against their
callable denominator rather than being dropped. Markers that are not callable
anywhere are named on the figure instead of being omitted silently.

Marker names and reference polarity come from stage 2's `*_marker_polarity.csv`,
so pfdhps A437G — where the 3D7 reference itself carries the mutant allele — is
reported as carriage, consistent with stage 3.

PNG (400 dpi), SVG and PDF are supplied with source CSVs.
""")
    print(f"Wrote six roadmap-driven publication figures to {a.outdir} (cohort n={n_total})")
if __name__=="__main__": main()

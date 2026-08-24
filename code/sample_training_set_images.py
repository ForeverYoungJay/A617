#!/usr/bin/env python3
"""Select representative 3D GB examples and render traceable QC panels."""
from __future__ import annotations
import argparse, csv, math
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from validate_3d_grain_graph import reconstruct

SETS = ("primary", "sensitivity", "review")

def rows(path):
    with path.open(newline="") as f: return list(csv.DictReader(f))

def number(x, key):
    try: return float(x[key])
    except (KeyError, ValueError): return math.nan

def choose(edges, targets, n=6):
    joined=[]
    t={(int(x["grain3d_a"]),int(x["grain3d_b"])):x for x in targets}
    for e in edges:
        key=(int(e["grain3d_a"]),int(e["grain3d_b"])); joined.append((key,e,t[key]))
    specs=[("largest_area",lambda x:number(x[1],"interface_area_um2"),True),
           ("smallest_area",lambda x:number(x[1],"interface_area_um2"),False),
           ("highest_misorientation",lambda x:number(x[1],"misorientation_deg"),True),
           ("shortest_lifetime",lambda x:number(x[1],"min_lifetime"),False)]
    finite=[x for x in joined if math.isfinite(number(x[2],"logratio_o_r3"))]
    if finite: specs[:0]=[("highest_O_r3",lambda x:number(x[2],"logratio_o_r3"),True),
                         ("lowest_O_r3",lambda x:number(x[2],"logratio_o_r3"),False)]
    else: specs += [("largest_volume_ratio",lambda x:number(x[1],"abs_log_volume_ratio"),True),
                    ("lowest_misorientation",lambda x:number(x[1],"misorientation_deg"),False)]
    out=[]; used=set()
    for reason,key,reverse in specs:
        candidates=[x for x in joined if x[0] not in used and math.isfinite(key(x))]
        if not candidates: continue
        item=sorted(candidates,key=key,reverse=reverse)[0]; used.add(item[0]); out.append((reason,*item))
        if len(out)==n: break
    for item in sorted(joined,key=lambda x:number(x[1],"interface_area_um2"),reverse=True):
        if item[0] not in used: out.append(("fallback_large_area",*item)); used.add(item[0])
    return out

def main():
    p=argparse.ArgumentParser(); p.add_argument("--train",type=Path,required=True); p.add_argument("--raw",type=Path,required=True); p.add_argument("--processed",type=Path,required=True); p.add_argument("--map",type=Path,required=True); p.add_argument("--out",type=Path,required=True); p.add_argument("--per-set",type=int,default=6); a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    section=defaultdict(lambda:defaultdict(list))
    for x in rows(a.map): section[int(x["grain3d_id"])][int(x["slice"])].append(int(x["local_grain_id"]))
    cache={}; manifest=[]
    for dataset in SETS:
        selected=choose(rows(a.train/dataset/"edges.csv"),rows(a.train/dataset/"targets.csv"),a.per_set)
        figures=[]
        for reason,(ga,gb),edge,target in selected:
            if len(figures)==a.per_set: break
            common=sorted(set(section[ga])&set(section[gb])); mid=len(common)//2; order=sorted(common,key=lambda z:abs(common.index(z)-mid))
            found=None
            for z in order:
                if z not in cache:
                    ang=next(a.raw.glob(f"*sliceimage-{z:03d}*.ang")); cache[z]=reconstruct(ang)
                rec=cache[z]; pairs=[tuple(sorted((x,y))) for x in section[ga][z] for y in section[gb][z]]
                pair=next((x for x in pairs if x in rec[4]),None)
                if pair is not None: found=(z,pair,rec); break
            if found is None: continue
            z,pair,(lab,q,valid,nodes,gb_edges,direct,metrics)=found; nr,nc=lab.shape
            ids=np.unique(gb_edges[pair][1]); rr,cc=np.divmod(ids,nc); pad=35
            r0,r1=max(0,int(rr.min())-pad),min(nr,int(rr.max())+pad+1); c0,c1=max(0,int(cc.min())-pad),min(nc,int(cc.max())+pad+1)
            rgb=np.abs(q[...,1:]); rgb/=np.maximum(rgb.sum(axis=2,keepdims=True),1e-12); rgb=np.where(valid[...,None],rgb,0)
            feature=a.processed/f"ebsd-sliceimage-{z:03d}_features.csv"; header=feature.open().readline().strip().split(","); ci_i=header.index("ci"); o_i=header.index("eds_ok_counts")
            ci=np.loadtxt(feature,delimiter=",",skiprows=1,usecols=ci_i).reshape(nr,nc); oxygen=np.loadtxt(feature,delimiter=",",skiprows=1,usecols=o_i).reshape(nr,nc); oxygen=np.log1p(np.maximum(oxygen,0)); lo,hi=np.nanpercentile(oxygen,[1,99])
            fig,axs=plt.subplots(1,3,figsize=(10.2,3.5),constrained_layout=True)
            images=(rgb,ci,oxygen); cmaps=(None,"gray","magma"); titles=("IPF-like","CI","O EDS (log)")
            for ax,img,cmap,title in zip(axs,images,cmaps,titles):
                ax.imshow(img[r0:r1,c0:c1],cmap=cmap,vmin=(lo if title.startswith("O") else None),vmax=(hi if title.startswith("O") else None),origin="upper")
                ax.scatter(cc-c0,rr-r0,s=3,c="#00ffff",linewidths=0); ax.set_title(title); ax.axis("off")
            fig.suptitle(f"{dataset.title()} · {reason} · GB {ga}–{gb} · slice {z}",fontsize=11)
            path=a.out/f"{dataset}_{reason}_gb{ga}-{gb}_slice{z:03d}.png"; fig.savefig(path,dpi=220); plt.close(fig); figures.append(path)
            manifest.append({"dataset":dataset,"selection_reason":reason,"grain3d_a":ga,"grain3d_b":gb,"slice":z,"local_grain_a":pair[0],"local_grain_b":pair[1],"edge_direct_in_slice":int(pair in direct),"interface_area_um2":edge["interface_area_um2"],"misorientation_deg":edge["misorientation_deg"],"min_lifetime":edge["min_lifetime"],"logratio_o_r3":target["logratio_o_r3"],"target_qc":target["target_qc"],"image":path.name})
        if figures:
            fig,axs=plt.subplots(2,3,figsize=(15,8));
            for ax,path in zip(axs.flat,figures): ax.imshow(plt.imread(path)); ax.axis("off")
            for ax in axs.flat[len(figures):]: ax.axis("off")
            fig.suptitle(f"{dataset.title()} representative QC samples",fontsize=15); fig.tight_layout(); fig.savefig(a.out/f"{dataset}_contact_sheet.png",dpi=180); plt.close(fig)
    with (a.out/"sample_manifest.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)
    assert all(sum(x["dataset"]==s for x in manifest)==a.per_set for s in SETS)
    print({s:sum(x["dataset"]==s for x in manifest) for s in SETS})

if __name__=="__main__": main()

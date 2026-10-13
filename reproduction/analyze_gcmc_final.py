"""Analyze existing final geometry with explicitly preset H depth and Cu descriptors."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1];initial=read(root/'GCMC/GCMC/POSCAR');final=read(a.directory/'CONTCAR')
initcu=initial[[i for i,x in enumerate(initial) if x.symbol=='Cu']];cu=final[[i for i,x in enumerate(final) if x.symbol=='Cu']]
assert len(cu)==len(initcu)
top_ids=np.flatnonzero(np.isclose(initcu.positions[:,2],initcu.positions[:,2].max(),atol=.01));fixed=float(np.median(initcu.positions[top_ids,2]));moving=float(np.median(cu.positions[top_ids,2]));h=final[[i for i,x in enumerate(final) if x.symbol=='H']]
r={'presets':{'top_layer_ids':'initial highest Cu layer, z tolerance 0.01 A','ambiguous_H_band_A':0.2,'Cu_cutoff_A':3.0,'GCN_denominator':12},'limitation':'H below a reference plane is a geometric diagnostic, not identified tetrahedral/octahedral Hsub. Cu layers can mix; no equilibrium or author exact GCN claim.','composition':dict(Counter(final.get_chemical_symbols())),'H_depth':{}}
for name,z in [('initial_top_plane',fixed),('moving_initial_top_ids',moving)]:
 dz=h.positions[:,2]-z;r['H_depth'][name]={'reference_z_A':z,'above_plus_0.2_A':int(sum(dz>.2)),'below_minus_0.2_A':int(sum(dz<-.2)),'within_band':int(sum(abs(dz)<=.2)),'H_z_minus_plane_A':dz.tolist()}
i,j=neighbor_list('ij',cu,3.0);cn=np.bincount(i,minlength=len(cu));gcn=np.zeros(len(cu));np.add.at(gcn,i,cn[j]/12.0)
r['Cu_coordination_histogram']=dict(Counter(map(int,cn)));r['top_Cu_GCN_mean']=float(gcn[top_ids].mean());r['top_Cu_CN_mean']=float(cn[top_ids].mean());r['top_Cu_z_displacement_range_A']=[float((cu.positions[top_ids,2]-initcu.positions[top_ids,2]).min()),float((cu.positions[top_ids,2]-initcu.positions[top_ids,2]).max())]
a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

"""Test height hypotheses on existing generated sites, without rerunning PH-SA."""
import json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.geometry import find_mic
root=Path(__file__).resolve().parents[1];initial=read(root/'GCMC/GCMC/POSCAR');top=float(initial.positions[:,2].max())
raw=json.loads((root/'reproduction/results/phsa_cu100/raw_sites.json').read_text());old=json.loads((root/'GCMC/GCMC/adsorption_sites.json').read_text());report={'scope':'Algebraic height transformation of existing surface sites; hypothesis, not author provenance. Internal sites and topology not regenerated.','initial_top_z_A':top,'species':{}}
for species,key,rnew,rtarget in [('H','H_surface',1.793,1.63),('CO','CO_surface',2.178,2.08)]:
 b=np.array([s['position'] for s in raw[key]]);h=b[:,2]-top;rad=h*h-rnew*rnew+rtarget*rtarget;fallback=rad<0;b[:,2]=top+np.where(fallback,h,np.sqrt(np.maximum(rad,0)))
 errors=[];normaldots=[]
 for i,point in enumerate(old[species]['positions']):
  _,d=find_mic(b-point,initial.cell,pbc=[True,True,False]);j=int(d.argmin());errors.append(float(d[j]))
  if species=='CO':normaldots.append(float(np.dot(old[species]['normals'][i],raw[key][j]['normal'])))
 report['species'][species]={'generated_effective_bond_A':rnew,'hypothesized_effective_bond_A':rtarget,'unchanged_fallback_sites':int(sum(fallback)),'max_nearest_distance_A':max(errors),'median_nearest_distance_A':float(np.median(errors)),'matched_within_0.0001_A':sum(x<.0001 for x in errors),'CO_min_normal_dot':min(normaldots) if normaldots else None}
(root/'reproduction/results/static_analysis/height_hypothesis.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

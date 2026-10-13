"""Static analysis only: existing sites and structures; no energy evaluation."""
import json
from pathlib import Path
from collections import Counter
import numpy as np
from ase.io import read
from ase.geometry import find_mic
from ase.neighborlist import neighbor_list
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'reproduction/results/static_analysis';out.mkdir(exist_ok=True)
atoms=read(ROOT/'GCMC/GCMC/POSCAR');cell=atoms.cell
raw=json.loads((ROOT/'reproduction/results/phsa_cu100/raw_sites.json').read_text())
old=json.loads((ROOT/'GCMC/GCMC/adsorption_sites.json').read_text())
report={}
for species,key in [('H','H_surface'),('CO','CO_surface')]:
 a=np.asarray(old[species]['positions']);b=np.asarray([s['position'] for s in raw[key]])
 rows=[]
 for i,p in enumerate(a):
  v=b-p;v[:,2]=0;_,d=find_mic(v,cell,pbc=[True,True,False]);j=int(d.argmin())
  row={'original_index':i,'generated_index':j,'xy_distance_A':float(d[j]),'delta_z_A':float(b[j,2]-p[2])}
  if species=='CO':row['normal_dot']=float(np.dot(old[species]['normals'][i],raw[key][j]['normal']))
  rows.append(row)
 report[species]={'xy_matches_0.01_A':sum(r['xy_distance_A']<.01 for r in rows),'unique_matched_generated':len(set(r['generated_index'] for r in rows)),'delta_z_groups_A':dict(Counter(round(r['delta_z_A'],5) for r in rows)),'rows':rows}
(out/'site_alignment.json').write_text(json.dumps(report,indent=2)+'\n')
paths=[ROOT/'GCMC/GCMC/POSCAR',ROOT/'reproduction/results/gcmc_baseline50_mu0/CONTCAR']+list((ROOT/'DFT calculation file').rglob('*'))
results=[]
for p in paths:
 if not p.is_file() or (p.name not in ['POSCAR','CONTCAR']) or ('GCMC simulation results' not in str(p) and p not in paths[:2]):continue
 a=read(p);cu=a[[i for i,x in enumerate(a) if x.symbol=='Cu']]
 i,j=neighbor_list('ij',cu,3.0);cn=np.bincount(i,minlength=len(cu));gcn=np.zeros(len(cu));np.add.at(gcn,i,cn[j]/12.0)
 results.append({'path':str(p.relative_to(ROOT)),'composition':dict(Counter(a.get_chemical_symbols())),'Cu_coordination_histogram':dict(Counter(map(int,cn))),'Cu_GCN_mean':float(gcn.mean()),'Cu_GCN_min':float(gcn.min()),'Cu_z_range_A':[float(cu.positions[:,2].min()),float(cu.positions[:,2].max())]})
(out/'structure_statistics.json').write_text(json.dumps({'presets':{'Cu_neighbor_cutoff_A':3.0,'GCN_bulk_normalization':12},'limitation':'Preset structural descriptors; not author exact GCN or equilibrium coverage.','structures':results},indent=2)+'\n')
print({k:{x:y for x,y in v.items() if x!='rows'} for k,v in report.items()});print('structures',len(results))

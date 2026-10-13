"""Summarize existing proposal and optimizer logs; never evaluate a potential."""
import argparse,json
from pathlib import Path
from collections import Counter
p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
def lines(name):
 f=a.directory/name
 if not f.exists():return []
 result=[]
 for line in f.read_text().splitlines():
  try:result.append(json.loads(line))
  except json.JSONDecodeError:break
 return result
m=lines('mc_steps.jsonl');o=lines('optimization_log.jsonl');s={'proposals':len(m),'accepted':sum(x['accepted'] for x in m),'operations':dict(Counter(x['operation'] for x in m)),'optimizations':len(o),'converged':sum(x['converged'] for x in o),'nonconverged_steps':[x['step'] for x in o if not x['converged']],'total_bfgs_steps':sum(x['bfgs_steps'] for x in o),'longest_bfgs':max((x['bfgs_steps'] for x in o),default=0),'over20':sum(x['bfgs_steps']>20 for x in o),'max_final_force_eV_per_A':max((x['max_force_eV_per_angstrom'] for x in o),default=None),'final':m[-1] if m else None,'output_file_bytes':sum(x.stat().st_size for x in a.directory.rglob('*') if x.is_file()),'blocks':[],'limitation':'Block means are diagnostics, not proof of convergence. Approximate gas chemical potentials and supplied acceptance rule.'}
for i in range(0,len(m),50):
 b=m[i:i+50];s['blocks'].append({'start_step':i,'count':len(b),'accepted':sum(x['accepted'] for x in b),'mean_n_CO':sum(x['n_CO'] for x in b)/len(b),'mean_n_H':sum(x['n_H'] for x in b)/len(b),'mean_energy_eV':sum(x['current_energy_eV'] for x in b)/len(b)})
if a.output:a.output.write_text(json.dumps(s,indent=2)+'\n')
print(json.dumps(s,indent=2))

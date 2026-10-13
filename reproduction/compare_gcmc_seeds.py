"""Compare existing GCMC records across seeds; no new proposals or model evaluation."""
import argparse,json
from pathlib import Path
from collections import Counter
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--first',type=Path,required=True);p.add_argument('--second',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();a.output_dir.mkdir(exist_ok=True)
series=[];report={'scope':'Same input/model/approximate mu; different random seeds. Early-chain diagnostics, not equilibrium estimates.','trajectories':[]}
for directory in [a.first,a.second]:
 rows=[json.loads(x) for x in (directory/'mc_steps.jsonl').read_text().splitlines()];config=json.loads((directory/'run_config.json').read_text());series.append((config['seed'],rows));blocks=[]
 for i in range(0,len(rows),50):
  b=rows[i:i+50];blocks.append({'start_step':i,'count':len(b),'mean_n_CO':sum(x['n_CO'] for x in b)/len(b),'mean_n_H':sum(x['n_H'] for x in b)/len(b),'accepted':sum(x['accepted'] for x in b)})
 operations={op:{'proposals':sum(x['operation']==op for x in rows),'accepted':sum(x['operation']==op and x['accepted'] for x in rows)} for op in sorted(set(x['operation'] for x in rows))}
 report['trajectories'].append({'seed':config['seed'],'proposals':len(rows),'accepted':sum(x['accepted'] for x in rows),'final_n_CO':rows[-1]['n_CO'],'final_n_H':rows[-1]['n_H'],'blocks':blocks,'operations':operations})
for key in ['input_sha256','mu_co','mu_h','temperature','opt_steps']:
 configs=[json.loads((d/'run_config.json').read_text()) for d in [a.first,a.second]]
 if configs[0][key]!=configs[1][key]:raise ValueError('Input mismatch: '+key)
(a.output_dir/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True)
for seed,rows in series:
 for ax,key in zip(axes,['n_CO','n_H']):ax.step([x['step']+1 for x in rows],[x[key] for x in rows],where='post',label=str(seed))
axes[0].set_ylabel('CO count');axes[1].set_ylabel('Total H count');axes[1].set_xlabel('Completed MC proposals')
for ax in axes:ax.legend(title='Seed');ax.grid(alpha=.25)
fig.suptitle('Cu(100): approximate chemical potentials; early-chain comparison');fig.tight_layout();fig.savefig(a.output_dir/'composition_vs_proposals.png',dpi=180);plt.close(fig)
print(json.dumps(report,indent=2))

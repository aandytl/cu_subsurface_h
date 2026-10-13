"""Scientific integration check: 4 continuous vs 2+2 resumed proposals. Requires permission."""
import argparse,json,os,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-root',required=True,type=Path);a=p.parse_args();a.output_root.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[1];env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','TF_INTRA_OP_PARALLELISM_THREADS']:env[k]='4'
env['TF_INTER_OP_PARALLELISM_THREADS']='1';env['CUDA_VISIBLE_DEVICES']='-1'
base=[sys.executable,str(root/'GCMC/GCMC/gcmc.py'),str(root/'GCMC/GCMC/POSCAR'),'--dp_model',str(root/'DP-GEN_file/DP-GEN_file/model_4.pb'),'--sites-file',str(root/'GCMC/GCMC/adsorption_sites.json'),'--mu-h','-2.4292938222817155','--mu-co','-14.28483274990397','--opt-steps','400','--seed','20261010']
for name,steps,resume in [('continuous',4,False),('split',2,False),('split',4,True)]:
 with (a.output_root/(name+str(steps)+'.log')).open('w') as log:
  subprocess.run(base+['--steps',str(steps),'--output-dir',str((a.output_root/name).resolve())]+(['--resume'] if resume else []),env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
x=json.loads((a.output_root/'continuous/results.json').read_text());y=json.loads((a.output_root/'split/results.json').read_text());equal=x==y
(a.output_root/'verification.json').write_text(json.dumps({'proposal_records_equal':equal,'scope':'4 continuous vs 2+2 resumed, same seed/input, 8 proposals total'},indent=2)+'\n')
if not equal:raise RuntimeError('Resume diverged from continuous execution')
print('Resume proposal records equal')

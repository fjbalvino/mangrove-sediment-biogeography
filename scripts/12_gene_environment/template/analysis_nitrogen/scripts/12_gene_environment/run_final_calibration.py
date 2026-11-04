"""Safe entry point for final calibration, including checkpoint provenance.

Never reuse simulation chunks after changing inputs, model code or the plan.
Use a fresh extracted package for a new functional matrix, and audit its panel
and sample design before changing the frozen reference assertions.
"""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);a=p.parse_args()
root=a.root.resolve();out=root/'calibration_final';scripts=root/'scripts/12_gene_environment'
inputs=['results/study_coordinate_groups.tsv','design_review/unit_design_manifest.tsv','results/group_detection_ge1.tsv','results/group_relative_representation.tsv','inputs/metadata_functional.csv']
code=['12G_PREPARE_FINAL_CALIBRATION.py','12G_RUN_FINAL_CALIBRATION.R']
paths=[root/x for x in inputs]+[scripts/x for x in code]
fingerprint={str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in paths}
guard=out/'resume_guard.json'
has_chunks=any((out/'chunks').glob('chunk_*.rds'))
if has_chunks and (not guard.exists() or json.loads(guard.read_text())!=fingerprint):
    sys.exit('Existing checkpoints lack matching provenance. Start in a fresh copy; do not mix changed inputs/code with old simulations.')
out.mkdir(exist_ok=True)
guard.write_text(json.dumps(fingerprint,indent=2)+'\n')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
subprocess.run([sys.executable,str(scripts/'12G_PREPARE_FINAL_CALIBRATION.py')],check=True,env=env)
subprocess.run(['Rscript',str(scripts/'12G_RUN_FINAL_CALIBRATION.R'),str(root),'validate'],check=True,env=env)
subprocess.run(['Rscript',str(scripts/'12G_RUN_FINAL_CALIBRATION.R'),str(root),'final'],check=True,env=env)
subprocess.run([sys.executable,str(scripts/'12H_SUMMARIZE_FINAL_CALIBRATION.py')],check=True,env=env)

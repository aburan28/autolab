#!/usr/bin/env sage -python
import argparse, datetime, hashlib, json
from pathlib import Path

def main():
    p=argparse.ArgumentParser(); p.add_argument('--primary',type=Path,required=True); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    x=json.loads(a.primary.read_text(encoding='ascii')); r=x['records']; g=x['gates']
    checks={'primary_preflight_pass':x['preflight_pass'],'finite_grid_retained':r['row_count']==324,'all_transition_gates_retained':all(g.values()),'mutation_control_retained':r['one_term_correction_mutation_changed_count']==324,'global_boundary_retained':'global cartier extension' in x['claim_status'].lower()}
    o={'schema':'ecdlp.h018.poincare-q-origin-transition.n609m.structural-replay.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'primary_sha256':hashlib.sha256(a.primary.read_bytes()).hexdigest(),'checks':checks,'pass':all(checks.values())}
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
    if not o['pass']: raise RuntimeError('N609M replay failed')
    print(json.dumps({'checks':sum(checks.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()

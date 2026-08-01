#!/usr/bin/env sage -python
"""N609Q: Q-series of the determinant's P-origin leading coefficient."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import vector
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture, source_corrections
from poincare_origin_boundary_n609c import constant_coefficient

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--out',type=Path,required=True); args=parser.parse_args()
    target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); corr=source_corrections(ring,u); determinant=matrix.det(); rows=[]
    for level in (0,1,17):
        moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9))
        regular7=moving[7]+moving[8]*corr[7]
        coefficient=-determinant*regular7
        rows.append({'level':level,'coefficient_valuation':int(coefficient.valuation()),'coefficient_constant':int(constant_coefficient(coefficient,field)),'regular7_valuation':int(regular7.valuation()),'regular7_constant':int(constant_coefficient(regular7,field))})
    gates={'three_levels_used':len(rows)==3,'determinant_is_q_origin_unit':determinant.valuation()==0,'p_origin_leading_coefficient_is_q_regular':all(row['coefficient_valuation']>=0 for row in rows),'p_origin_leading_coefficient_is_nonzero_at_crossing':all(row['coefficient_constant']!=0 for row in rows),'regular7_boundary_is_uniform_45':all(row['regular7_constant']==45 for row in rows)}
    o={'schema':'ecdlp.h018.poincare-crossing-coefficient.n609q.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / DETERMINANT P-ORIGIN LEADING COEFFICIENT EXTENDS THROUGH Q-ORIGIN / MODEL-BOUND / TOY-EVIDENCE / FULL_TWO_VARIABLE_CECH_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'formula':'coefficient of v^8 D_t in corrected Q frame is -det(V_Q)*r_7(Q)'},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The determinant P-origin leading coefficient is Q-origin regular and nonzero for t=0,1,17; its corrected section coefficient has uniform Q-boundary value 45. This gives formal compatibility of the determinant P-origin trivialization with the Q-origin frame at the rational crossing.','next_requirement':'Construct the full two-variable overlap, including the corrected Cauchy basis and non-rational charts, before any global Cartier, smoothness, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claim.'}
    args.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
    if not o['preflight_pass']: raise RuntimeError('N609Q crossing coefficient failed')
    print(json.dumps({'checks':sum(gates.values()),'output':str(args.out)},sort_keys=True))
if __name__=='__main__': main()

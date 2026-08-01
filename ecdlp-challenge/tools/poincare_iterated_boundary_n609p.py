#!/usr/bin/env sage -python
"""N609P: iterated Q-origin then P-origin H018 boundary check."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import GF
from poincare_global_sections_n608r import P
from poincare_origin_boundary_n609c import constant_coefficient, regular_sections

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--out',type=Path,required=True); args=parser.parse_args()
    _target, field, sections=regular_sections(); base=GF(P); rows=[]
    for level in (0,1,17):
        coeff=[base(constant_coefficient(sections['x'][i]-base(level)*sections['1'][i],field)) for i in range(9)]
        # In the L(9O) P-origin frame v^9, x^3 y is the constant term and x^4
        # is the v coefficient.  N609C's x^3y coefficient must vanish.
        rows.append({'level':level,'q_boundary_x3y_coefficient':int(coeff[8]),'q_boundary_x4_coefficient':int(coeff[7]),'iterated_p_origin_constant':int(coeff[8]),'iterated_p_origin_linear_coefficient':int(coeff[7]),'simple_p_origin_zero':coeff[8]==0 and coeff[7]!=0})
    gates={'three_selected_levels_used':len(rows)==3,'q_origin_transition_removes_x3y_term':all(row['q_boundary_x3y_coefficient']==0 for row in rows),'iterated_boundary_has_uniform_nonzero_linear_term':all(row['iterated_p_origin_linear_coefficient']==45 for row in rows),'all_members_have_simple_iterated_boundary_zero':all(row['simple_p_origin_zero'] for row in rows)}
    o={'schema':'ecdlp.h018.poincare-iterated-boundary.n609p.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / ITERATED Q-ORIGIN TO P-ORIGIN H018 BOUNDARY COMPATIBILITY / MODEL-BOUND / TOY-EVIDENCE / TWO_VARIABLE_CECH_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'frame_order':'N609M corrected Q-origin frame, then the L(9O) P-origin frame v^9'},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'After the corrected Q-origin boundary, every selected member has zero constant and uniform nonzero v-linear coefficient 45 in the P-origin L(9O) frame. Thus the iterated boundary has a simple rational crossing zero.','next_requirement':'Construct a genuine two-variable overlap and non-rational chart transition; iterated boundary data alone does not prove a global Cech cocycle, Cartier divisor, smoothness, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claim.'}
    args.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
    if not o['preflight_pass']: raise RuntimeError('N609P iterated-boundary preflight failed')
    print(json.dumps({'checks':sum(gates.values()),'output':str(args.out)},sort_keys=True))
if __name__=='__main__': main()

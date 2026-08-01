#!/usr/bin/env sage -python
"""N609R: audit the corrected-basis contribution to the P-origin lead."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import vector
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture, source_corrections
from poincare_origin_boundary_n609c import constant_coefficient

def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    target,cover,pi,phi=fixture(); field,ring,u,_support,matrix,xvalues=origin_data(target,cover,pi,phi); corr=source_corrections(ring,u); det=matrix.det(); rows=[]
    # B_8^reg=u^-8(b_8-sum_i corr_i b_i) has x^4 coefficient -u^-8*corr_7=-u^-1.
    for level in (0,1,17):
        moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); reg7=moving[7]+moving[8]*corr[7]; reg8=moving[8]*u**8
        corrected_x4=reg7-reg8/u
        raw_x4=moving[7]
        rows.append({'level':level,'basis_x4_correction_identity':corrected_x4==raw_x4,'full_leading_valuation':int((-det*corrected_x4).valuation()),'full_leading_constant':int(constant_coefficient(-det*corrected_x4,field)),'shortcut_leading_constant':int(constant_coefficient(-det*reg7,field))})
    gates={'three_levels_used':len(rows)==3,'corrected_cauchy_x4_term_cancels_reg7_correction':all(row['basis_x4_correction_identity'] for row in rows),'full_leading_is_q_regular':all(row['full_leading_valuation']>=0 for row in rows),'full_leading_is_nonzero_at_crossing':all(row['full_leading_constant']!=0 for row in rows)}
    o={'schema':'ecdlp.h018.poincare-crossing-full-leading-audit.n609r.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / FULL CORRECTED-BASIS P-ORIGIN LEADING COEFFICIENT AUDIT / MODEL-BOUND / TOY-EVIDENCE / TWO_VARIABLE_CECH_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'rows':rows,'formula':'full x^4 coefficient is r7-r8/u=moving7'},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'After including the corrected Cauchy basis, its x^4 contribution exactly cancels the r8 correction to r7. The full P-origin leading coefficient is -det(V_Q)*moving_7; its Q-regularity and crossing value are recorded directly.','next_requirement':'Use the audited full coefficient, not the r7 shortcut, in any two-variable Cech construction; extend through remaining basis orders and non-rational charts before global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}
    a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii')
    if not o['preflight_pass']: raise RuntimeError('N609R full leading audit failed')
    print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()

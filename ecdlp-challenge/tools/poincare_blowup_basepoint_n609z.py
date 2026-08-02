#!/usr/bin/env sage -python
"""N609Z: test whether the direct exceptional root is a common pencil base point."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import vector
from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture
def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(); p.add_argument('--out',type=Path,required=True); a=p.parse_args(); target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); formal=target.formal_group(); slope=field(83); vals={}
 for level in (0,1,17):
  moving=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); parameter=ring(slope)*u; point=(ring(formal.x(32)(parameter)),ring(formal.y(32)(parameter))); raw=sum(moving[i]*moving_basis(point,support)[i] for i in range(9)); vals[level]=-det*parameter**8*raw
 rows=[]
 for level in (0,1,17):
  diff=vals[level]-vals[0]; rows.append({'level':level,'normalized_valuation':int(vals[level].valuation()),'level_difference_valuation':None if diff.is_zero() else int(diff.valuation())})
 gates={'all_levels_have_order_four_contact':all(row['normalized_valuation']==4 for row in rows),'all_nonzero_level_differences_vanish_beyond_contact_order':all(row['level']==0 or row['level_difference_valuation'] is None or row['level_difference_valuation']>4 for row in rows)}
 o={'schema':'ecdlp.h018.poincare-blowup-basepoint.n609z.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / H018 DIRECT EXCEPTIONAL ROOT IS A COMMON INFINITELY-NEAR PENCIL BASEPOINT CANDIDATE / MODEL-BOUND / TOY-EVIDENCE / BASEPOINT_RESOLUTION_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM','records':{'slope':'83','rows':rows},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'The direct exceptional root has common order-four contact across the selected pencil, and the first level differences have the recorded higher order. This is local evidence for an infinitely-near common base point, not yet a resolved global base locus.','next_requirement':'Compute the complete local ideal and blow-up multiplicity sequence, then remove any proven common base contribution before genus, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}; a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii');
 if not o['preflight_pass']: raise RuntimeError('N609Z basepoint probe failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()

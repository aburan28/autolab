#!/usr/bin/env sage -python
"""N610A: local base-ideal leading terms at the H018 exceptional point."""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
from sage.all import vector
from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture
def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(); p.add_argument('--out',type=Path,required=True); a=p.parse_args(); target,cover,pi,phi=fixture(); field,ring,u,support,matrix,xvalues=origin_data(target,cover,pi,phi); det=matrix.det(); formal=target.formal_group(); slope=field(83); values={}
 for level in (0,1,17):
  m=matrix.solve_right(xvalues-ring(level)*vector(ring,[1]*9)); q=ring(slope)*u; point=(ring(formal.x(32)(q)),ring(formal.y(32)(q))); values[level]=-det*q**8*sum(m[i]*moving_basis(point,support)[i] for i in range(9))
 rows=[]
 for level in (1,17):
  d=values[level]-values[0]; rows.append({'level':level,'difference_valuation':int(d.valuation()),'difference_leading_coefficient':str(d[d.valuation()]),'order_eight_unit':d.valuation()==8 and d[8]!=0})
 gates={'both_nonzero_level_differences_are_order_eight_units':all(r['order_eight_unit'] for r in rows),'base_member_retains_order_four_contact':values[0].valuation()==4}
 o={'schema':'ecdlp.h018.poincare-blowup-baseideal.n610a.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'OBSERVATION / LOCAL H018 EXCEPTIONAL PENCIL BASE-IDEAL LEADING TERMS / MODEL-BOUND / TOY-EVIDENCE / MULTIPLICITY_SEQUENCE_AND_GLOBAL_RESOLUTION_OPEN / NO_ECDLP_CLAIM','records':{'slope':'83','rows':rows,'base_member_valuation':int(values[0].valuation())},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'At the exceptional point, the selected pencil has one member with order-four normal contact and both independent level differences are order-eight units. This supports the local base-ideal candidate (x+O(u^4),u^8 unit) once the slope coordinate derivative is included.','next_requirement':'Compute the slope derivative and successive blow-up transforms of this local ideal, then subtract only a proven base scheme before using genus, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.'}; a.out.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n',encoding='ascii');
 if not o['preflight_pass']: raise RuntimeError('N610A baseideal probe failed')
 print(json.dumps({'checks':sum(gates.values()),'output':str(a.out)},sort_keys=True))
if __name__=='__main__': main()

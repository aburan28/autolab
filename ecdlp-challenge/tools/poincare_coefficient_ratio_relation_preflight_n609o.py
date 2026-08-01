#!/usr/bin/env sage -python
"""N609O: nonlinear coefficient-ratio target-relation screen for H018."""
from __future__ import annotations

import argparse, datetime, hashlib, json
from pathlib import Path

from sage.all import GF, EllipticCurve

from poincare_linear_projection_n608v import P, coefficient_table, complete_sources


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--out',type=Path,required=True); args=parser.parse_args()
    field=GF(P); curve=EllipticCurve(field,[0,0,0,1,24]); coefficients=coefficient_table(curve); levels=(0,1,17)
    sources={level:complete_sources(curve,coefficients,level) for level in levels}
    rows=[]
    for denominator in (7,8):
        if not all(vector[denominator] for vector in coefficients.values()):
            continue
        for numerator in range(9):
            scalars={q:int(vector[numerator]/vector[denominator]) for q,vector in coefficients.items()}
            image_sizes=[]
            for level in levels:
                values={point+scalars[base]*base for point,base in sources[level]}
                image_sizes.append(len(values))
            rows.append({'numerator':numerator,'denominator':denominator,'minimum_level_image_size':min(image_sizes),'maximum_level_image_size':max(image_sizes),'singleton_level_count':sum(size==1 for size in image_sizes),'mean_level_image_size':sum(image_sizes)/len(image_sizes)})
    best=min(rows,key=lambda row:(row['maximum_level_image_size'],row['mean_level_image_size'],row['numerator'],row['denominator']))
    gates={'all_admitted_ratio_maps_screened':len(rows)==18,'all_three_certified_levels_used':len(sources)==3,'no_ratio_level_map_is_single_valued':all(row['singleton_level_count']==0 for row in rows),'best_ratio_map_retains_nontrivial_level_ambiguity':best['maximum_level_image_size']>1}
    output={'schema':'ecdlp.h018.poincare-coefficient-ratio-relation-preflight.n609o.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'NEGATIVE RESULT / H018 COEFFICIENT-RATIO MAPS DO NOT GIVE A SINGLE-VALUED TARGET RELATION / MODEL-BOUND / TOY-EVIDENCE / OTHER_NONLINEAR_RELATIONS_OPEN / NO_ECDLP_CLAIM','records':{'certified_levels':list(levels),'map_rows':rows,'best_ratio_map':best},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'For every admitted raw moving-frame coefficient ratio a(Q)=c_i(Q)/c_7(Q) or c_i(Q)/c_8(Q), none of the certified complete H018 levels determines P+a(Q)Q. This compact nonlinear correction family supplies no direct target-bearing relation law.','next_requirement':'Test a nonlinear correction carrying a proved composition law or derive one from a compact normalization/Jacobian before factor-base, rank, descent, cost, or ECDLP claims.'}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',encoding='ascii')
    if not output['preflight_pass']: raise RuntimeError('N609O coefficient-ratio preflight failed')
    print(json.dumps({'checks':sum(gates.values()),'output':str(args.out)},sort_keys=True))
if __name__=='__main__': main()

#!/usr/bin/env sage -python
"""N609N: additive target-bearing relation screen for complete H018 levels."""
from __future__ import annotations

import argparse, datetime, hashlib, json
from pathlib import Path

from sage.all import GF, EllipticCurve

from poincare_linear_projection_n608v import P, coefficient_table, complete_sources


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--out',type=Path,required=True); args=parser.parse_args()
    field=GF(P); curve=EllipticCurve(field,[0,0,0,1,24]); coefficients=coefficient_table(curve)
    levels=(0,1,17)
    level_sources={level:complete_sources(curve,coefficients,level) for level in levels}
    rows=[]
    for multiplier in range(int(curve.cardinality())):
        image_sizes=[]
        for level, sources in level_sources.items():
            values={point + multiplier*base for point,base in sources}
            image_sizes.append(len(values))
        rows.append({'multiplier':multiplier,'minimum_level_image_size':min(image_sizes),'maximum_level_image_size':max(image_sizes),'singleton_level_count':sum(size==1 for size in image_sizes),'mean_level_image_size':sum(image_sizes)/len(image_sizes)})
    best=min(rows,key=lambda row:(row['maximum_level_image_size'],row['mean_level_image_size'],row['multiplier']))
    gates={'all_109_additive_maps_screened':len(rows)==109,'all_three_certified_levels_used':len(level_sources)==3,'no_additive_level_map_is_single_valued':all(row['singleton_level_count']==0 for row in rows),'best_map_retains_nontrivial_level_ambiguity':best['maximum_level_image_size']>1}
    output={'schema':'ecdlp.h018.poincare-additive-relation-preflight.n609n.v1','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha256_file(Path(__file__)),'claim_status':'NEGATIVE RESULT / H018 CERTIFIED COMPLETE LEVELS DO NOT GIVE A SINGLE-VALUED ADDITIVE TARGET RELATION / MODEL-BOUND / TOY-EVIDENCE / NONLINEAR_RELATIONS_OPEN / NO_ECDLP_CLAIM','records':{'curve_order':int(curve.cardinality()),'certified_levels':list(levels),'rows':rows,'best_additive_map':best},'gates':gates,'preflight_pass':all(gates.values()),'strongest_valid_statement':'For every c in E(F_103), none of the three certified complete H018 level correspondences has P+cQ uniquely determined across all declared sources. Thus the raw complete-level data supplies no direct additive target-bearing relation P+cQ=H(t) on the tested levels.','next_requirement':'Test a genuinely nonlinear correspondence with an explicit compact correction state, or prove a relation law from a compact normalization/Jacobian before factor-base, rank, descent, cost, or ECDLP claims.'}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',encoding='ascii')
    if not output['preflight_pass']: raise RuntimeError('N609N additive relation preflight failed')
    print(json.dumps({'checks':sum(gates.values()),'output':str(args.out)},sort_keys=True))
if __name__=='__main__': main()

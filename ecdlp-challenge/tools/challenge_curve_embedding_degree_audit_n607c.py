#!/usr/bin/env sage -python
"""N607C: exact MOV embedding-degree audit for the 40-bit public curve."""
import argparse,json
from pathlib import Path
from sage.all import Mod,factor
P=616883774851;N=616882790773
def main():
 p=argparse.ArgumentParser();p.add_argument("--output",type=Path,required=True);a=p.parse_args();k=int(Mod(P,N).multiplicative_order());fs=str(factor(N-1));checks={"factorization":fs=="2^2 * 3^2 * 17135633077","order_divides_n_minus_one":(N-1)%k==0,"exact_large_embedding_degree":k==17135633077};o={"schema":"ecdlp.challenge-curve.embedding-degree-audit.n607c.v1","claim_status":"RESTRICTED THEOREM / EXACT_MOV_EMBEDDING_DEGREE / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM","records":{"p":P,"n":N,"n_minus_one_factorization":fs,"embedding_degree":k},"gates":checks,"preflight_pass":all(checks.values()),"next_requirement":"Any pairing route must name a different target subgroup or additional special structure; this MOV embedding degree is not practical."};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"checks":sum(checks.values()),"output":str(a.output)}))
if __name__=="__main__":main()

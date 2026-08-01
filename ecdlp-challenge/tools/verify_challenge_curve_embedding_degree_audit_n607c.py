#!/usr/bin/env sage -python
import argparse,json
from pathlib import Path
from sage.all import Mod
def main():
 p=argparse.ArgumentParser();p.add_argument("--primary",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=json.loads(a.primary.read_text());k=int(Mod(616883774851,616882790773).multiplicative_order());c={"primary":d["preflight_pass"] is True,"order":k==17135633077,"boundary":k>200};o={"verified":all(c.values()),"checks":c};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"verified":o["verified"],"output":str(a.output)}))
if __name__=="__main__":main()

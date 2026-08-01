#!/usr/bin/env sage -python
import argparse,json
from pathlib import Path
def main():
 p=argparse.ArgumentParser();p.add_argument("--primary",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();d=json.loads(a.primary.read_text());c={"primary":d["preflight_pass"] is True,"ranks":[x["rank"] for x in d["rows"]]==[9,9,9],"boundary":d["global_h018_sections_constructed"] is False};o={"verified":all(c.values()),"checks":c};a.output.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n");print(json.dumps({"verified":o["verified"],"output":str(a.output)}))
if __name__=="__main__":main()

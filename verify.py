#!/usr/bin/env python3
"""Independently replay published coordinate and exclusion certificates."""
import argparse,concurrent.futures,hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def replay(item):
 label,checker,certificate,timeout=item
 try:
  p=subprocess.run([sys.executable,str(checker),str(certificate)],capture_output=True,text=True,timeout=timeout)
  report=json.loads(p.stdout)
  if checker.parent.name=='bi-certifier':
   passed=p.returncode==0 and report.get('status')=='EXCLUDED_NONFLAT' and all(report.get(k) is True for k in ('root_certified','convex_metric_certified','excluded'))
  else:
   passed=p.returncode==0 and report.get('verified') is True
  return dict(item=label,passed=passed,report=report)
 except Exception as error:return dict(item=label,passed=False,error=str(error))

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--body',type=int,choices=range(1,199),metavar='1..198')
 parser.add_argument('--jobs',type=int,default=1)
 parser.add_argument('--timeout',type=float,default=900)
 args=parser.parse_args()
 if args.jobs<1 or args.timeout<=0:parser.error('jobs and timeout must be positive')
 catalogue=json.loads((ROOT/'catalogue.json').read_text());items=[]
 for row in catalogue:
  n=row['body_id']
  if args.body is not None and args.body!=n:continue
  certificate=ROOT/'certificates'/f'{n:03d}'/'coordinates.json'
  if hashlib.sha256(certificate.read_bytes()).hexdigest()!=row['sha256']:
   print(f'Body {n}: certificate checksum mismatch');return 2
  items.append((f'Body {n}',ROOT/'software/metric-converter/check.py',certificate,args.timeout))
 if args.body is None:
  for certificate in sorted((ROOT/'exclusions').glob('*/certificate.json')):
   items.append((f'Exclusion {certificate.parent.name}',ROOT/'software/bi-certifier/check.py',certificate,args.timeout))
 failed=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
  for result in pool.map(replay,items):
   print(json.dumps(result),flush=True)
   if not result['passed']:failed.append(result['item'])
 print(json.dumps({'checked':len(items),'passed':len(items)-len(failed),'failed':failed}))
 return 2 if failed else 0
if __name__=='__main__':raise SystemExit(main())

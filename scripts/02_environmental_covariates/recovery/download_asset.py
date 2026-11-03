from pathlib import Path
import urllib.request,json,sys,hashlib,time
u=sys.argv[1];p=Path(sys.argv[2]);p.parent.mkdir(exist_ok=True,parents=True)
if p.exists():print('cached',p,p.stat().st_size,flush=True);sys.exit()
t=time.time();h=hashlib.sha256()
with urllib.request.urlopen(u,timeout=45) as r,p.with_suffix(p.suffix+'.partial').open('wb') as f:
 headers=dict(r.headers);n=0
 while c:=r.read(2**20):
  f.write(c);h.update(c);n+=len(c)
  if n%(100*2**20)==0:print(p.name,n,'bytes',round(time.time()-t),flush=True)
 if 'Content-Length' in headers:assert n==int(headers['Content-Length'])
p.with_suffix(p.suffix+'.partial').rename(p)
p.with_suffix('.download.json').write_text(json.dumps({'url':u,'headers':headers,'sha256':h.hexdigest(),'bytes':n,'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())},indent=2))
print('done',p,n,round(time.time()-t),flush=True)

import pathlib,re,urllib.request,concurrent.futures,hashlib,csv,datetime
root=pathlib.Path('.')
jobs=[]
for acc in ['GSE263995','GSE289825']:
 s=(root/'metadata/single_cell'/f'{acc}_family.soft').read_text(encoding='utf-8')
 if acc=='GSE263995':
  wanted=['GSM8208510','GSM8208511','GSM8391107','GSM8391110','GSM8391113','GSM8391116']
  for block in s.split('^SAMPLE = ')[1:]:
   if block.splitlines()[0] in wanted:
    jobs.extend((acc,u.replace('ftp://','https://')) for u in re.findall(r'!Sample_supplementary_file_\d+ = (\S+)',block))
 else: jobs.extend((acc,u.replace('ftp://','https://')) for u in re.findall(r'!Series_supplementary_file = (\S+)',s))
def get(job):
 acc,url=job; f=root/'data/raw/single_cell'/acc/url.split('/')[-1]
 if not f.exists():
  with urllib.request.urlopen(url,timeout=120) as response, open(str(f)+'.part','wb') as out:
   while b:=response.read(1024*1024):out.write(b)
  pathlib.Path(str(f)+'.part').rename(f)
 return dict(GSM=f.name.split('_')[0] if f.name.startswith('GSM') else '',accession=acc,url=url,local_relative_path=f.as_posix(),download_date=datetime.datetime.now(datetime.timezone.utc).isoformat(),file_size_bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
rows=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for row in pool.map(get,jobs):
  rows.append(row); print(row['local_relative_path'],row['file_size_bytes'],flush=True)
with open(root/'metadata/single_cell/download_manifest.tsv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


"""Audited sector-boundary sensitivity."""
import argparse
import csv
import json
import math
import statistics
import struct
import zipfile
from pathlib import Path
from reproduce_beta import records,cell,SECTORS

ROOT=Path(__file__).resolve().parents[1]
YEARS=list(range(2000,2015))
CODES=[c for s in ('R','e','m','x','s') for c in SECTORS[s]]
FINAL=('CONS_h','CONS_np','CONS_g','GFCF')

def finite(v,where):
    if not isinstance(v,(int,float)) or not math.isfinite(v):raise ValueError(f'Missing/non-numeric {where}')
    return float(v)

def read_wiot(path):
    with zipfile.ZipFile(path) as z:
        ss=z.read('xl/sharedStrings.bin');strings=[]
        for kind,p,end in records(ss):
            if kind==19:
                n=struct.unpack_from('<I',ss,p+1)[0]
                strings.append(ss[p+5:p+5+2*n].decode('utf-16le'))
        data=z.read('xl/worksheets/sheet1.bin')
    hi={};hc={};buyers={};final={};sources={};chfinal={};va=None
    row=-1;code='';country='';values={};fvalues={}
    def flush():
        nonlocal va
        if code=='VA':
            if va is not None:raise ValueError('Duplicate VA')
            va=values.copy()
        if code in SECTORS['e']+SECTORS['R']:
            key=(code,country)
            if key in sources:raise ValueError('Duplicate source row')
            sources[key]=values.copy()
        if code in CODES and country=='CHN':
            if code in chfinal:raise ValueError('Duplicate China industry')
            chfinal[code]=fvalues.copy()
    for kind,p,end in records(data):
        if kind==0:
            if row>=6:flush()
            row=struct.unpack_from('<I',data,p)[0];code='';country='';values={};fvalues={}
            if row==5:
                buyers={c:hi[c] for c in hi if hc.get(c)=='CHN' and hi[c] in CODES}
                final={c:hi[c] for c in hi if hc.get(c)=='CHN' and hi[c] in FINAL}
                if sorted(buyers.values())!=sorted(CODES) or sorted(final.values())!=sorted(FINAL):raise ValueError('China header coverage mismatch')
            continue
        if kind not in range(1,12):continue
        col=struct.unpack_from('<I',data,p)[0]
        if row in (2,4) and col>=4:(hi if row==2 else hc)[col]=str(cell(data,kind,p,strings)).strip()
        elif row>=6:
            if col==0:code=str(cell(data,kind,p,strings)).strip()
            elif col==2:country=str(cell(data,kind,p,strings)).strip()
            elif col in buyers and (code in SECTORS['e']+SECTORS['R'] or code=='VA'):
                values[col]=finite(cell(data,kind,p,strings),(row,col))
            elif col in final and code in CODES and country=='CHN':
                fvalues[final[col]]=finite(cell(data,kind,p,strings),(row,col))
    flush()
    countries={k[1] for k in sources}
    if len(countries)!=44 or 'CHN' not in countries or set(sources)!={(c,ct) for c in SECTORS['e']+SECTORS['R'] for ct in countries}:raise ValueError('Global supplier coverage mismatch')
    if va is None or set(va)!=set(buyers) or set(chfinal)!=set(CODES):raise ValueError('VA/final-use coverage mismatch')
    if any(set(v)!=set(buyers) for v in sources.values()) or any(set(v)!=set(FINAL) for v in chfinal.values()):raise ValueError('Selected cell coverage mismatch')
    result=[]
    for col,c in buyers.items():
        r={'code':c,'VA':va[col]}
        for name,codes in [('processing',SECTORS['e']),('raw',SECTORS['R'])]:
            r[name+'_global']=math.fsum(v[col] for (sc,ct),v in sources.items() if sc in codes)
        r.update(chfinal[c]);result.append(r)
    return result

def read_sea(path):
    import openpyxl
    w=openpyxl.load_workbook(path,read_only=True,data_only=True)
    rows=w['DATA'].iter_rows(values_only=True);h=next(rows);ix={str(v):i for i,v in enumerate(h)}
    if not all(k in ix for k in ['country','variable','code']+[str(y) for y in YEARS]):raise ValueError('SEA headers mismatch')
    found={}
    for r in rows:
        if r[ix['country']]!='CHN' or r[ix['code']] not in CODES or r[ix['variable']] not in ('CAP','COMP'):continue
        key=(r[ix['code']],r[ix['variable']])
        if key in found:raise ValueError('Duplicate SEA key')
        found[key]=[finite(r[ix[str(y)]],(key,y)) for y in YEARS]
    w.close()
    if set(found)!={(c,v) for c in CODES for v in ('CAP','COMP')}:raise ValueError('SEA industry coverage mismatch')
    return [dict(year=y,code=c,CAP=found[(c,'CAP')][i],COMP=found[(c,'COMP')][i]) for i,y in enumerate(YEARS) for c in CODES]

def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load_csv(path):
    with path.open(newline='') as f:
        return [{k:int(v) if k=='year' else v if k=='code' else float(v) for k,v in r.items()} for r in csv.DictReader(f)]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--wiot-dir',type=Path);ap.add_argument('--sea',type=Path)
    ap.add_argument('--data-dir',type=Path,default=ROOT/'data/measurement')
    ap.add_argument('--output-dir',type=Path)
    ap.add_argument('--mapping-file',type=Path)
    a=ap.parse_args();out=a.output_dir or a.data_dir;out.mkdir(parents=True,exist_ok=True)
    wp=a.data_dir/'sector_wiot_annual.csv';sp=a.data_dir/'sector_sea_annual.csv'
    if a.wiot_dir:
        wi=[]
        for y in YEARS:
            wi.extend(dict(year=y,**r) for r in read_wiot(a.wiot_dir/f'WIOT{y}_Nov16_ROW.xlsb'))
        write_csv(out/'sector_wiot_annual.csv',wi)
    else:wi=load_csv(wp)
    if a.sea:
        sea=read_sea(a.sea);write_csv(out/'sector_sea_annual.csv',sea)
    else:sea=load_csv(sp)
    def keyed(rows):
        d={(r['year'],r['code']):r for r in rows}
        if len(d)!=len(rows) or set(d)!={(y,c) for y in YEARS for c in CODES}:raise ValueError('Need exactly 15 years x 53 industries')
        return d
    wd=keyed(wi);sd=keyed(sea)
    choices=json.loads((a.mapping_file or a.data_dir/'sector_mapping_choices.json').read_text())
    ids=[c['id'] for c in choices['cases']]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate case identifiers')
    diagnostic=[]
    for c in SECTORS['x']:
        cons=math.fsum(wd[(y,c)]['CONS_h']+wd[(y,c)]['CONS_np']+wd[(y,c)]['CONS_g'] for y in YEARS)
        inv=math.fsum(wd[(y,c)]['GFCF'] for y in YEARS)
        yearly=[]
        for y in YEARS:
            r=wd[(y,c)];C=math.fsum(r[k] for k in FINAL if k!='GFCF');I=r['GFCF']
            if min(C,I)<0:raise ValueError('Negative selected final use')
            if C+I>0:yearly.append(I/(C+I))
        diagnostic.append(dict(code=c,consumption=cons,GFCF=inv,GFCF_share=inv/(cons+inv) if cons+inv>0 else None,
          mean_annual_GFCF_share=statistics.mean(yearly) if yearly else None,
          consumption_dominant_years=sum(v<0.5 for v in yearly),positive_final_use_years=len(yearly)))
    results={};annual=[]
    for case in choices['cases']:
        mapping=case['mapping']
        flat=[c for s in ('R','e','m','x','s') for c in mapping[s]]
        if len(flat)!=len(set(flat)) or set(flat)!=set(CODES):raise ValueError('Mapping must partition original 53 industries')
        if mapping['R']!=SECTORS['R'] or mapping['e']!=SECTORS['e']:raise ValueError('This reader supports fixed raw/processing supplier mapping only')
        shares=[];lab={}
        for s in ('m','s','x','e'):
            group=mapping[s];betas=[];cap=comp=0.
            for y in YEARS:
                VA=math.fsum(wd[(y,c)]['VA'] for c in group)
                inp=math.fsum(wd[(y,c)]['raw_global' if s=='e' else 'processing_global'] for c in group)
                CAP=math.fsum(sd[(y,c)]['CAP'] for c in group);COMP=math.fsum(sd[(y,c)]['COMP'] for c in group)
                if min(VA,inp,CAP,COMP)<0 or min(VA+inp,CAP+COMP)<=0:raise ValueError('Invalid sector aggregates')
                beta=inp/(VA+inp);betas.append(beta);cap+=CAP;comp+=COMP;lab[(s,y)]=COMP
                annual.append(dict(case_id=case['id'],year=y,sector=s,VA=VA,input_global=inp,CAP=CAP,COMP=COMP,beta=beta))
            b=statistics.mean(betas);alpha=(1-b)*cap/(cap+comp)
            if min(alpha,b)<=0 or alpha+b>=1:raise ValueError('Infeasible production shares')
            shares.append([alpha,b])
        s0=lab[('s',2000)]/(lab[('s',2000)]+lab[('m',2000)])
        s14=lab[('s',2014)]/(lab[('s',2014)]+lab[('m',2014)])
        results[case['id']]={'shares':shares,'targets':[1.63,s0,s14-s0],'description':case['description'],'mapping':mapping}
    report={'years':YEARS,'sector_order':['m','s','x','e'],'aggregation':'Original pooled CAP/(CAP+COMP), arithmetic mean of annual global beta; full precision','final_use_scope':'CHN producer rows to CHN final-use columns; excludes inventories, exports and intermediate uses','final_use_diagnostic':diagnostic,'cases':results}
    (out/'sector_mapping_parameters.json').write_text(json.dumps(report,indent=2))
    write_csv(out/'sector_mapping_annual.csv',annual)
    lines=['function d = sector_mapping_parameter_data()','% Generated by build_sector_mapping_parameters.py; edit mapping JSON then regenerate.']
    for name,r in results.items():
        if not name.isidentifier() or not name.isascii():raise ValueError('Invalid case identifier')
        lines.append('d.'+name+'.shares=['+'; '.join(' '.join(format(v,'.17g') for v in row) for row in r['shares'])+'];')
        lines.append('d.'+name+'.targets=['+'; '.join(format(v,'.17g') for v in r['targets'])+'];')
    lines.append('end');(out/'sector_mapping_parameter_data.m').write_text('\n'.join(lines)+'\n')
    print('Mapping inputs saved:', out)

if __name__=='__main__':main()

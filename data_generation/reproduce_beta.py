"""Read selected WIOT XLSB cells."""
import argparse
import csv
import json
import math
import statistics
import struct
import zipfile
from pathlib import Path

SECTORS = {
 'R':['B'], 'e':['C19','C20','C23','C24','D35'],
 'm':['C10-C12','C13-C15','C16','C17','C18','C21','C22','C31_C32'],
 'x':['C25','C26','C27','C28','C29','C30','C33','F'],
 's':['E36','E37-E39','G45','G46','G47','H49','H50','H51','H52','H53',
 'I','J58','J59_J60','J61','J62_J63','K64','K65','K66','L68','M69_M70',
 'M71','M72','M73','M74_M75','N','O84','P85','Q','R_S','T','U']}

def records(data):
    p=0
    while p<len(data):
        kind=0; shift=0
        while True:
            c=data[p];p+=1;kind|=(c&127)<<shift;shift+=7
            if not c&128:break
        size=0;shift=0
        while True:
            c=data[p];p+=1;size|=(c&127)<<shift;shift+=7
            if not c&128:break
        end=p+size
        if end>len(data):raise ValueError('Truncated BIFF12 record')
        yield kind,p,end
        p=end

def cell(data,kind,p,strings):
    if kind==7:return strings[struct.unpack_from('<I',data,p+8)[0]]
    if kind in (5,9):return struct.unpack_from('<d',data,p+8)[0]
    if kind==2:
        raw=struct.unpack_from('<I',data,p+8)[0]
        if raw&2:
            signed=struct.unpack('<i',struct.pack('<I',raw))[0]
            value=signed>>2
        else:value=struct.unpack('<d',struct.pack('<II',0,raw&0xfffffffc))[0]
        return value/100 if raw&1 else value
    if kind in (6,8):
        n=struct.unpack_from('<I',data,p+8)[0]
        return data[p+12:p+12+2*n].decode('utf-16le')
    if kind==1:return None
    raise ValueError(f'Unsupported selected cell record {kind}')

def read_year(path):
    with zipfile.ZipFile(path) as z:
        ss=z.read('xl/sharedStrings.bin')
        strings=[]
        for kind,p,end in records(ss):
            if kind==19:
                n=struct.unpack_from('<I',ss,p+1)[0]
                strings.append(ss[p+5:p+5+2*n].decode('utf-16le'))
        data=z.read('xl/worksheets/sheet1.bin')
    header_ind={};header_cty={};buyers={};source_rows=[];va=None
    current=-1;rowcode='';rowcty='';values={};retain=False
    selected_suppliers=set(SECTORS['e']+SECTORS['R'])
    def flush():
        nonlocal va
        if rowcode=='VA':
            if va is not None:raise ValueError('Duplicate VA row')
            va=values.copy()
        elif rowcode in selected_suppliers:source_rows.append((rowcode,rowcty,values.copy()))
    for kind,p,end in records(data):
        if kind==0:
            if current>=6:flush()
            current=struct.unpack_from('<I',data,p)[0]
            rowcode='';rowcty='';values={};retain=False
            if current==5:
                for s in ('m','x','s','e'):
                    cols={c:header_ind[c] for c in header_ind if header_cty.get(c)=='CHN' and header_ind[c] in SECTORS[s]}
                    if sorted(cols.values())!=sorted(SECTORS[s]):raise ValueError(f'Buyer coverage mismatch {s}')
                    buyers[s]=list(cols)
                wanted=set(c for cols in buyers.values() for c in cols)
            continue
        if kind not in range(1,12):continue
        col=struct.unpack_from('<I',data,p)[0]
        if current in (2,4) and col>=4:
            value=cell(data,kind,p,strings)
            (header_ind if current==2 else header_cty)[col]=str(value).strip()
        elif current>=6:
            if col==0:
                rowcode=str(cell(data,kind,p,strings)).strip()
                retain=rowcode in selected_suppliers or rowcode=='VA'
            elif col==2 and retain:rowcty=str(cell(data,kind,p,strings)).strip()
            elif retain and col in wanted:
                value=cell(data,kind,p,strings)
                if value is None or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError(f'Missing/non-numeric {current},{col}')
                values[col]=value
    flush()
    if va is None or set(va)!=wanted:raise ValueError('VA coverage mismatch')
    countries={cty for code,cty,vals in source_rows}
    if len(countries)!=44 or 'CHN' not in countries:raise ValueError('Supplier country coverage mismatch')
    if len(source_rows)!=44*6:raise ValueError('Supplier row count mismatch')
    keys={(code,cty) for code,cty,vals in source_rows}
    if len(keys)!=len(source_rows):raise ValueError('Duplicate supplier')
    for code,cty,vals in source_rows:
        if set(vals)!=wanted:raise ValueError(f'Incomplete selected supplier row {code} {cty}')
    out=[]
    for s in ('m','x','s','e'):
        cols=buyers[s];supplier=set(SECTORS['R'] if s=='e' else SECTORS['e'])
        VA=math.fsum(va[c] for c in cols)
        domestic=math.fsum(vals[c] for code,cty,vals in source_rows if code in supplier and cty=='CHN' for c in cols)
        global_=math.fsum(vals[c] for code,cty,vals in source_rows if code in supplier for c in cols)
        if min(VA,domestic,global_)<0 or VA+domestic<=0:raise ValueError('Invalid aggregate denominator/flow')
        out.append(dict(sector=s,VA=VA,input_domestic=domestic,input_import=global_-domestic,input_global=global_,beta_domestic=domestic/(VA+domestic),beta_global=global_/(VA+global_)))
    return out,dict(buyer_columns=len(wanted),supplier_rows=len(source_rows),supplier_countries=len(countries),selected_numeric_cells=(len(source_rows)+1)*len(wanted))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
    rows=[];checks={}
    for year in range(2000,2015):
        path=args.input_dir/f'WIOT{year}_Nov16_ROW.xlsb'
        if not path.is_file():raise FileNotFoundError(path)
        results,check=read_year(path);checks[str(year)]=check
        rows.extend(dict(year=year,**r) for r in results)
    summary={}
    config={'m':.3024,'x':.5415,'s':.1035,'e':.4737}
    for s in config:
        group=[r for r in rows if r['sector']==s]
        summary[s]={k:statistics.mean(r[k] for r in group) for k in ('beta_domestic','beta_global')}
        summary[s]['config_beta']=config[s];summary[s]['global_round4_matches_config']=round(summary[s]['beta_global'],4)==config[s]
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/'beta_annual.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    (args.output_dir/'beta_audit.json').write_text(json.dumps(dict(years=list(range(2000,2015)),checks=checks,summary=summary),indent=2))
    print('Input shares saved:', args.output_dir)

if __name__=='__main__':main()

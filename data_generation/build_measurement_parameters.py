"""Build editable MATLAB production-share inputs for measurement robustness."""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ORDER=('m','s','x','e')
SECTORS={
 'e':['C19','C20','C23','C24','D35'],
 'm':['C10-C12','C13-C15','C16','C17','C18','C21','C22','C31_C32'],
 'x':['C25','C26','C27','C28','C29','C30','C33','F'],
 's':['E36','E37-E39','G45','G46','G47','H49','H50','H51','H52','H53',
 'I','J58','J59_J60','J61','J62_J63','K64','K65','K66','L68','M69_M70',
 'M71','M72','M73','M74_M75','N','O84','P85','Q','R_S','T','U']}

def read_sea(path):
    import openpyxl
    w=openpyxl.load_workbook(path,read_only=True,data_only=True)
    it=w['DATA'].iter_rows(values_only=True);head=next(it)
    ix={str(v):j for j,v in enumerate(head)}
    required=['country','variable','code']+[str(y) for y in range(2000,2015)]
    if not all(k in ix for k in required):raise ValueError('SEA header mismatch')
    found={}
    for row in it:
        if row[ix['country']]!='CHN' or row[ix['variable']] not in ('CAP','COMP','LAB','VA'):continue
        key=(row[ix['variable']],row[ix['code']])
        if key in found:raise ValueError(f'Duplicate SEA row {key}')
        vals=[row[ix[str(y)]] for y in range(2000,2015)]
        if any(v is None or not isinstance(v,(int,float)) or not math.isfinite(v) for v in vals):raise ValueError(f'Missing SEA observation {key}')
        found[key]=vals
    out=[]
    for s in ORDER:
        for y in range(2000,2015):
            r=dict(year=y,sector=s)
            for v in ('CAP','COMP','LAB','VA'):
                r[v]=math.fsum(found[(v,c)][y-2000] for c in SECTORS[s])
            out.append(r)
    w.close();return out

def load_csv(path):
    with path.open(newline='') as f:rows=list(csv.DictReader(f))
    return [dict((k, int(v) if k=='year' else v if k=='sector' else float(v)) for k,v in r.items()) for r in rows]

def check_keys(rows):
    keys=[(r['year'],r['sector']) for r in rows]
    if len(set(keys))!=len(keys) or set(keys)!={(y,s) for y in range(2000,2015) for s in ORDER}:raise ValueError('Need exactly 15 years x 4 sectors')

def matrix_literal(matrix):return '['+'; '.join(' '.join(format(v,'.17g') for v in row) for row in matrix)+']'

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sea',type=Path)
    ap.add_argument('--wiot-dir',type=Path)
    ap.add_argument('--annual-beta',type=Path,default=ROOT/'data/measurement/beta_annual.csv')
    ap.add_argument('--annual-capital',type=Path,default=ROOT/'data/measurement/sea_capital_annual.csv')
    ap.add_argument('--output-dir',type=Path,default=ROOT/'data/measurement')
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    capital=read_sea(a.sea) if a.sea else load_csv(a.annual_capital)
    if a.wiot_dir:
        from reproduce_beta import read_year
        beta=[]
        for y in range(2000,2015):
            rows,_=read_year(a.wiot_dir/f'WIOT{y}_Nov16_ROW.xlsb')
            beta.extend(dict(year=y,**r) for r in rows)
    else:beta=load_csv(a.annual_beta)
    check_keys(capital);check_keys(beta)
    c={(r['year'],r['sector']):r for r in capital};b={(r['year'],r['sector']):r for r in beta}
    annual=[];summary={};rounded_beta=dict(zip(ORDER,[.3024,.1035,.5415,.4737]))
    for s in ORDER:
        cr=[c[(y,s)] for y in range(2000,2015)]
        pooled=math.fsum(r['CAP'] for r in cr)/math.fsum(r['CAP']+r['COMP'] for r in cr)
        vals=[]
        for y in range(2000,2015):
            r=c[(y,s)];v=b[(y,s)]
            if min(r['CAP'],r['COMP'])<=0:raise ValueError('Nonpositive factor income')
            tilde=r['CAP']/(r['CAP']+r['COMP'])
            for source in ('global','domestic'):
                inp=v['input_'+source];beta_value=inp/(v['VA']+inp)
                if abs(beta_value-v['beta_'+source])>1e-12:raise ValueError('Beta numerator/denominator mismatch')
            row=dict(year=y,sector=s,alpha_tilde=tilde,beta_global=v['beta_global'],
                     beta_domestic=v['beta_domestic'],alpha_global=(1-v['beta_global'])*tilde,
                     alpha_domestic=(1-v['beta_domestic'])*tilde)
            vals.append(row);annual.append(row)
        global_mean=statistics.mean(r['beta_global'] for r in vals)
        if round(global_mean,4)!=rounded_beta[s]:raise ValueError(f'Baseline beta reproduction failed {s}')
        domestic_mean=statistics.mean(r['beta_domestic'] for r in vals)
        summary[s]=dict(pooled_alpha_tilde=pooled,
            reproduced_alpha_using_rounded_beta=(1-rounded_beta[s])*pooled,
            annual_alpha=statistics.mean(r['alpha_global'] for r in vals),annual_beta=global_mean,
            domestic_alpha=(1-domestic_mean)*pooled,domestic_beta=domestic_mean)
    cases={
      'annual_average':[[summary[s]['annual_alpha'],summary[s]['annual_beta']] for s in ORDER],
      'domestic_inputs':[[summary[s]['domestic_alpha'],summary[s]['domestic_beta']] for s in ORDER]}
    for matrix in cases.values():
        if any(not (0<x[0] and 0<x[1] and sum(x)<1) for x in matrix):raise ValueError('Infeasible production shares')
    for name,rows in [('sea_capital_annual.csv',capital),('beta_annual.csv',beta),('production_shares_annual.csv',annual)]:
        with (a.output_dir/name).open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    provenance=dict(years=list(range(2000,2015)),sector_order=list(ORDER),sector_mapping=SECTORS,
      notes=['Annual: mean of each complete annual alpha/beta pair; global inputs.',
             'Domestic: original pooled capital-income aggregation; mean annual domestic beta.',
             'Original baseline remains the rounded values in solver_config.m.',
             'CAP and COMP in million local currency; WIOT VA/inputs in million USD. Only within-source ratios are combined.'],
      summary=summary,cases=cases)
    (a.output_dir/'measurement_parameters.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    text=['function d = measurement_parameter_data()',
      '% Generated by data_generation/build_measurement_parameters.py; edit scenarios in config.',
      "d.sector_order={'m','s','x','e'};",'d.years=2000:2014;']
    for key,matrix in cases.items():text.append(f'd.{key}={matrix_literal(matrix)};')
    text.append('end')
    (a.output_dir/'measurement_parameter_data.m').write_text('\n'.join(text)+'\n')
    print('Measurement inputs saved:', a.output_dir)

if __name__=='__main__':main()

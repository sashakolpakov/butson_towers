"""Freeze, launch and resume exact native Williamson axis searches.

Random sampling and annealing evaluate each axis word with its uniquely
forced signs. A positive result is rebuilt at every actual residue/shift.
The authoritative checkpoint retains the full MT19937-64 state and the
current annealing word; resumed workers do not replay the seed prefix.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import product
import json
from math import isqrt
import os
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def row_targets(data):
    """Complete finite exact DC/count targets within the supplied chart."""
    p,m,n=data['p'],data['m'],data['N'];result=set()
    if data['action']=='C3':
        for head in range(-n,n+1):
            for tail in range(-n,n+1):
                if abs(head)+abs(tail)>n or (n-head-tail)%2:continue
                a=1-m*head+m*tail;b=1+m*head+(m//3)*tail
                if a*a+3*b*b==4*p:result.add((a,b,b,b))
    elif data['action']=='identity':
        R=isqrt(4*p);choices=[r for r in range(-R,R+1) if (r-1)%m==0]
        for row in product(choices,repeat=4):
            if sum(r*r for r in row)!=4*p or (sum(row)-4)%(2*m):continue
            total=(sum(row)-4)//(2*m)
            if any((total-(r-1)//m)%2 for r in row):continue
            nets=tuple((total-(r-1)//m)//2 for r in row)
            if sum(nets)!=total:raise ArithmeticError('row conversion failed')
            if sum(map(abs,nets))<=n and (n-total)%2==0:result.add(row)
    else:raise ValueError('coordinate action must be identity or C3')
    if not result:raise ValueError('the chart has no compatible row norm/count target')
    return [list(r) for r in sorted(result)]


def serialize(data,targets):
    """No cached index arrays are needed: native startup rebuilds every orbit."""
    if data['N']!=(data['p']-1)//data['m'] or len(data['rows'])!=data['N']:
        raise ValueError('coefficient dimensions failed')
    lines=['DENSE_WILLIAMSON_AXES_V1',
           f"{data['p']} {data['m']} {data['mu']} {data['N']} {int(data['action']=='C3')}",str(len(targets))]
    lines += [' '.join(map(str,row)) for row in targets]
    for row in data['rows']:
        lines.append(str(len(row)))
        lines += [' '.join(map(str,term)) for term in row]
    return ('\n'.join(lines)+'\n').encode()


def prepare(args):
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    data=json.loads(args.coefficients.read_text());targets=row_targets(data)
    source=Path(__file__).resolve();cpp=source.with_suffix('.cpp')
    raw_input=serialize(data,targets)
    specification={'schema':'dense-williamson-native-specification-v1',
        'p':data['p'],'m':data['m'],'mu':data['mu'],'N':data['N'],'action':data['action'],
        'allowed_row_sums':targets,'mode':args.mode,'seed':args.seed,
        'maximum_witnesses':args.max_witnesses,'initial_temperature':float(args.temperature),
        'minimum_temperature':float(args.minimum_temperature),'epoch_steps':args.epoch_steps,
        'source_sha256':{'coefficients.json':sha(args.coefficients),
            'native_axis_search.cpp':sha(cpp),'native_axis_search.py':sha(source),
            'input.txt':hashlib.sha256(raw_input).hexdigest()},
        'checkpoint_contract':'Full mt19937_64 state, current/best axes, cumulative elapsed time and counters.',
        'exclusion_scope':'Only complete base-four enumeration is exhaustive. Random/annealing counters include repeated words.'}
    fingerprint=hashlib.sha256(json.dumps(specification,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    frozen=output/'frozen';parameters=output/'parameters.json'
    if parameters.exists():
        record=json.loads(parameters.read_text())
        if record['fingerprint']!=fingerprint or record['specification']!=specification:
            raise ValueError('resume requires identical immutable source/configuration identity')
        if sha(frozen/'source-manifest.json')!=record['frozen_manifest_sha256']:
            raise ValueError('frozen source manifest differs from its immutable parameter record')
        manifest=json.loads((frozen/'source-manifest.json').read_text())
        actual={str(p.relative_to(frozen)) for p in frozen.rglob('*') if p.is_file() and p.name!='source-manifest.json'}
        if set(manifest)!=actual:raise ValueError('frozen file set changed')
        for name,digest in manifest.items():
            path=frozen/name
            if path.is_symlink() or sha(path)!=digest:raise ValueError('frozen source hash mismatch: '+name)
    else:
        if frozen.exists():raise ValueError('partial frozen output exists; use a new directory')
        frozen.mkdir()
        for original,name in [(args.coefficients,'coefficients.json'),(cpp,'native_axis_search.cpp'),(source,'native_axis_search.py')]:
            shutil.copy2(original,frozen/name)
        (frozen/'input.txt').write_bytes(raw_input)
        compiler=shutil.which('c++')
        if compiler is None:raise RuntimeError('a C++17 compiler is required')
        subprocess.run([compiler,'-O3','-std=c++17',str(frozen/'native_axis_search.cpp'),'-o',str(frozen/'engine')],check=True)
        manifest={str(p.relative_to(frozen)):sha(p) for p in sorted(frozen.iterdir()) if p.is_file()}
        write_json(frozen/'source-manifest.json',manifest)
        write_json(parameters,{'fingerprint':fingerprint,'specification':specification,
            'created_utc':utc(),'initial_budget_seconds':args.seconds,
            'frozen_manifest_sha256':sha(frozen/'source-manifest.json')})
        for path in frozen.iterdir():path.chmod(0o555 if path.name=='engine' else 0o444)
        parameters.chmod(0o444)
    command=[str(frozen/'engine'),str(frozen/'input.txt'),str(output),str(args.seconds),str(args.seed),
        args.mode,fingerprint,str(args.max_patterns),str(args.max_witnesses),str(args.temperature),
        str(args.minimum_temperature),str(args.epoch_steps),str(args.nice),'1']
    resume=[sys.executable,'-B',str(frozen/'native_axis_search.py'),'--coefficients',str(frozen/'coefficients.json'),
        '--output',str(output),'--mode',args.mode,'--seed',str(args.seed),'--seconds',str(args.seconds),
        '--max-patterns',str(args.max_patterns),'--max-witnesses',str(args.max_witnesses),
        '--temperature',str(args.temperature),'--minimum-temperature',str(args.minimum_temperature),
        '--epoch-steps',str(args.epoch_steps),'--nice',str(args.nice)]
    return output,command,resume,fingerprint


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coefficients',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mode',choices=('random','anneal','enumerate'),default='random')
    parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--seed',type=int,default=20261127)
    parser.add_argument('--max-patterns',type=int,default=0,help='cumulative limit; zero means no pattern limit')
    parser.add_argument('--max-witnesses',type=int,default=1,help='zero keeps searching after guarded witnesses')
    parser.add_argument('--temperature',type=float,default=1000)
    parser.add_argument('--minimum-temperature',type=float,default=1)
    parser.add_argument('--epoch-steps',type=int,default=10000)
    parser.add_argument('--nice',type=int,default=0)
    parser.add_argument('--detach',action='store_true')
    args=parser.parse_args(argv)
    if (args.seconds<=0 or args.seed<0 or args.seed>=1<<64 or args.max_patterns<0 or args.max_witnesses<0
        or args.temperature<=0 or args.minimum_temperature<=0 or args.epoch_steps<=0 or not 0<=args.nice<=19):
        parser.error('invalid native resource/stream parameters')
    output,command,resume,fingerprint=prepare(args)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    if args.detach:
        with (output/'worker.log').open('ab',buffering=0) as log:
            process=subprocess.Popen(command,cwd=output,env=env,stdin=subprocess.DEVNULL,stdout=log,
                stderr=subprocess.STDOUT,start_new_session=True)
        write_json(output/'launch.json',{'pid':process.pid,'command':command,'resume_command':resume,
            'fingerprint':fingerprint,'launched_utc':utc(),'nice':args.nice})
        print(json.dumps({'pid':process.pid,'output':str(output),'fingerprint':fingerprint}))
    else:
        write_json(output/'launch.json',{'command':command,'resume_command':resume,
            'fingerprint':fingerprint,'launched_utc':utc(),'nice':args.nice})
        subprocess.run(command,cwd=output,env=env,check=True)


if __name__=='__main__':main()

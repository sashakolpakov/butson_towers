"""Native exact-count local search for twisted periodic T-matrix designs.

A nonzero objective is a search checkpoint, never an existence certificate.
Candidates with zero objective are expanded and verified at every shift.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
from experiments.scalar_lift import hadamard_frame_search as frame


def write_tables(path, data, head_net, tail_net):
    """Same exact coefficient model as the small engine, without word limits."""
    n, count = data['n'], len(data['rows'])
    h = np.zeros((n, 8, count), dtype=np.int16)
    k = np.zeros((n, n, 8, 8, count), dtype=np.int16)
    axes, signs = np.arange(8)//2, 1-2*(np.arange(8)%2)
    for ell, row in enumerate(data['rows']):
        for (i,j,delta), multiplicity in row.items():
            rotated = np.where(axes == 0, 0, 1+(axes-1+delta)%3)
            if i == j:
                h[i,:,ell] += multiplicity*(axes == rotated)
            else:
                values = multiplicity*(axes[:,None] == rotated[None,:])*signs[:,None]*signs[None,:]
                k[i,j,:,:,ell] += values
                k[j,i,:,:,ell] += values.T
    with Path(path).open('wb') as output:
        np.array([data['p'],data['m'],data['mu'],n,count,
                  1+data['m']*head_net,data['m']//3*tail_net],dtype=np.int32).tofile(output)
        h.tofile(output); k.tofile(output)
    return h,k


def verify_checkpoint(record, data, head_net, tail_net):
    states=record['states']
    b=frame.expand(states,data)
    sums=b.sum(axis=1).tolist()
    expected=[1+data['m']*head_net]+[data['m']//3*tail_net]*3
    if sums != expected:
        raise ValueError('checkpoint changed exact row sums')
    correlations=[int(np.sum(b*np.roll(b,-s,axis=1))) for s in data['shifts']]
    if correlations != record['correlations']:
        raise ValueError('compressed correlation table disagrees with expanded integers')
    if sum(c*c for c in correlations) != record['best_energy']:
        raise ValueError('checkpoint objective does not equal exact correlation energy')
    full=[int(np.sum(b*np.roll(b,-s,axis=1))) for s in range(data['p'])]
    zero=full == [data['p']]+[0]*(data['p']-1)
    if zero != record['exact_correlation_zero']:
        raise ValueError('zero-status assertion disagrees with all shifts')
    return {'all_expanded_correlations_checked':True,'exact_row_sums_checked':True,
            'T_witness_found':zero,'maximum_absolute_nonzero_shift':max(map(abs,full[1:]),default=0)}


def build(directory):
    compiler=shutil.which('c++') or shutil.which('clang++')
    if not compiler:raise RuntimeError('C++17 compiler needed')
    binary=Path(directory)/'t_sequence_local_search'
    snapshot=Path(directory)/'native-source.cpp'
    snapshot.write_bytes(Path(__file__).with_suffix('.cpp').read_bytes())
    subprocess.run([compiler,'-O3','-std=c++17',str(snapshot),'-o',str(binary)],check=True)
    return binary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--p',type=int,default=223);parser.add_argument('--m',type=int,default=3)
    parser.add_argument('--mu',type=int,default=183)
    parser.add_argument('--head-net',type=int,default=-5);parser.add_argument('--tail-net',type=int,default=3)
    parser.add_argument('--heads',type=int,default=5);parser.add_argument('--seconds',type=float,default=60)
    parser.add_argument('--seed',type=int,default=20261125223)
    parser.add_argument('--algorithm',choices=('tabu','anneal'),default='tabu')
    parser.add_argument('--free-heads',action='store_true',
                        help='allow balanced head/tail pair conversions while retaining exact coordinate sums')
    parser.add_argument('--absolute-weight',type=int,default=0)
    parser.add_argument('--parity-weight',type=int,default=0)
    parser.add_argument('--initial',type=Path)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--workdir',type=Path)
    args=parser.parse_args()
    if args.absolute_weight<0 or args.parity_weight<0:parser.error("objective weights must be nonnegative")
    if args.seconds<=0:parser.error("search duration must be positive")
    target=[1+args.m*args.head_net]+[args.m//3*args.tail_net]*3
    if sum(x*x for x in target)!=args.p:
        parser.error("target row sums do not satisfy the necessary squared-norm identity")
    folder=args.workdir or Path(tempfile.mkdtemp(prefix='t-local-'))
    folder.mkdir(parents=True,exist_ok=True);args.report.parent.mkdir(parents=True,exist_ok=True)
    data=frame.chart(args.p,args.m,args.mu)
    table=folder/'coefficients.bin';write_tables(table,data,args.head_net,args.tail_net)
    runner_snapshot=folder/'runner-source.py'
    runner_snapshot.write_bytes(Path(__file__).read_bytes())
    binary=build(folder)
    source_hash=hashlib.sha256((folder/'native-source.cpp').read_bytes()).hexdigest()
    subprocess.run([str(binary),str(table),str(args.report),str(args.seconds),str(args.seed),str(args.heads),'check',str(args.absolute_weight),str(args.parity_weight),'-',str(int(args.free_heads))],check=True)
    command=[str(binary),str(table),str(args.report),str(args.seconds),str(args.seed),str(args.heads),args.algorithm,str(args.absolute_weight),str(args.parity_weight)]
    if args.initial:
        initial=json.loads(args.initial.read_text())['states']
        initial_b=frame.expand(initial,data)
        expected=[1+args.m*args.head_net]+[args.m//3*args.tail_net]*3
        if initial_b.sum(axis=1).tolist()!=expected:
            raise ValueError('initial states have incompatible exact row sums')
        if sum(state//2==0 for state in initial[1:])!=args.heads:
            raise ValueError('initial states have incompatible head count')
        initial_file=folder/'initial-states.txt'
        initial_file.write_text(' '.join(map(str,initial))+'\n')
        command.append(str(initial_file))
    if not args.initial:command.append('-')
    command.append(str(int(args.free_heads)))
    subprocess.run(command,check=True)
    record=json.loads(args.report.read_text())
    record['independent_expanded_verification']=verify_checkpoint(record,data,args.head_net,args.tail_net)
    record.update(algorithm=args.algorithm,head_orbits=sum(s//2==0 for s in record['states'][1:]),
                  initial_head_orbits=args.heads,free_head_count=args.free_heads,
                  head_net=args.head_net,tail_net=args.tail_net,
                  source_sha256=source_hash,
                  runner_source_sha256=hashlib.sha256(runner_snapshot.read_bytes()).hexdigest(),
                  coefficient_table_sha256=hashlib.sha256(table.read_bytes()).hexdigest(),
                  native_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  completed_new_order_claimed=False)
    even_path=Path(str(args.report)+'.even.json')
    if even_path.exists():
        even=json.loads(even_path.read_text())
        even['independent_expanded_verification']=verify_checkpoint(even,data,args.head_net,args.tail_net)
        even.update(head_orbits=sum(s//2==0 for s in even['states'][1:]),
                    initial_head_orbits=args.heads,free_head_count=args.free_heads,
                    head_net=args.head_net,tail_net=args.tail_net,
                    source_sha256=source_hash,completed_new_order_claimed=False)
        even_path.write_text(json.dumps(even,indent=2)+'\n')
    args.report.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['independent_expanded_verification']))

if __name__=='__main__':main()

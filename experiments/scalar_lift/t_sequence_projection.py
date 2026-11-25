"""Douglas--Rachford search for exact twisted T-matrix designs.

The Fourier and assignment steps use floating point only as search guidance.
Every saved candidate is independently expanded using integers. Only an exact
zero at every nonzero shift is labelled a completed T-design.
"""
from __future__ import annotations
import argparse
import ctypes
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import numpy as np
from experiments.scalar_lift import hadamard_frame_search as base


def build_assignment(directory):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).with_suffix('.cpp')
    snapshot=directory/'assignment-source.cpp';snapshot.write_bytes(source.read_bytes())
    runner=directory/'projection-source.py';runner.write_bytes(Path(__file__).read_bytes())
    library=directory/'t_sequence_assignment.dylib'
    compiler=shutil.which('c++') or shutil.which('clang++')
    if compiler is None:raise RuntimeError('C++17 compiler is required')
    subprocess.run([compiler,'-O3','-std=c++17','-shared','-fPIC',str(snapshot),'-o',str(library)],check=True)
    lib=ctypes.CDLL(str(library.resolve()))
    function=lib.t_sequence_assign
    function.argtypes=[ctypes.c_int,np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS'),
                       np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS')]
    function.restype=ctypes.c_int
    free_function=lib.t_sequence_assign_free_heads
    free_function.argtypes=[ctypes.c_int,np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),
                            ctypes.c_int,ctypes.c_int,
                            np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS')]
    free_function.restype=ctypes.c_int
    function.free_heads=free_function
    return function,lib


class Projections:
    def __init__(self,data,heads,head_net,tail_net,assignment):
        self.data=data;self.assignment=assignment
        self.heads=heads;self.head_net=head_net;self.tail_net=tail_net
        n=data['n']-1
        if heads is None:
            if (n+head_net+tail_net)%2 or abs(head_net)+abs(tail_net)>n:
                raise ValueError('requested exact DC sums are impossible')
            self.counts=None
        else:
            tails=n-heads
            counts=[(heads+head_net)//2,(heads-head_net)//2,
                    (tails+tail_net)//2,(tails-tail_net)//2]
            if heads<abs(head_net) or tails<abs(tail_net) or (heads+head_net)%2 or (tails+tail_net)%2:
                raise ValueError('requested exact sign counts are impossible')
            self.counts=np.array(counts,dtype=np.int32)
            if min(counts)<0 or sum(counts)!=n:raise ValueError('incorrect category counts')
        self.target=np.array([1+data['m']*head_net]+[data['m']//3*tail_net]*3,dtype=np.int64)
        if int(self.target@self.target)!=data['p']:raise ValueError('target DC vector has incorrect norm')
        self.positions=np.array([[x for x,k in orbit] for orbit in data['orbits'][1:]],dtype=np.intp)
        phase=np.array([[k for x,k in orbit] for orbit in data['orbits'][1:]],dtype=np.intp)
        self.coordinates=np.stack([np.zeros_like(phase) if axis==0 else 1+(axis-1+phase)%3 for axis in range(4)])
        self.radius=float(np.sqrt(data['p']))

    def assign_scores(self,scores):
        scores=np.ascontiguousarray(scores,dtype=np.float64)
        n=self.data['n']-1
        if scores.shape!=(n,4):raise ValueError('numeric score dimensions differ')
        out=np.empty(n,dtype=np.int32)
        status=(self.assignment.free_heads(n,scores,self.head_net,self.tail_net,out)
                if self.counts is None else self.assignment(n,scores,self.counts,out))
        if status:raise ValueError(f'assignment failed with status{status}')
        return [0]+out.tolist()

    def discrete(self,values):
        values=np.asarray(values,dtype=np.float64)
        if values.shape!=(4,self.data['p']):raise ValueError('wrong iterate dimensions')
        scores=np.stack([values[self.coordinates[a],self.positions].sum(axis=1) for a in range(4)],axis=1)
        states=self.assign_scores(scores)
        return states,base.expand(states,self.data)

    def spectral(self,frequencies):
        frequencies=np.asarray(frequencies,dtype=np.complex128)
        if frequencies.shape!=(4,self.data['p']):raise ValueError('wrong spectral dimensions')
        norms=np.sqrt(np.sum(np.abs(frequencies)**2,axis=0))
        answer=frequencies*(self.radius/np.where(norms>1e-14,norms,1.0))
        zero=norms<=1e-14
        answer[:,zero]=0;answer[0,zero]=self.radius
        positive=np.arange(1,(self.data['p']+1)//2)
        answer[:,self.data['p']-positive]=answer[:,positive].conj()
        answer[:,0]=self.target
        return answer


def verify(states,data,target):
    b=base.expand(states,data)
    if not np.array_equal(b.sum(axis=1),target):raise ValueError('exact sign counts failed')
    full=[int(np.sum(b*np.roll(b,-s,axis=1))) for s in range(data['p'])]
    classes=[full[s] for s in data['shifts']]
    return {'states':states,'correlations':classes,'row_sums':b.sum(axis=1).tolist(),
            'best_energy':sum(c*c for c in classes),
            'maximum_absolute_nonzero_shift':max(map(abs,full[1:]),default=0),
            'all_expanded_shift_correlations_checked':True,
            'exact_T_design':full==[data['p']]+[0]*(data['p']-1)}


def controls(assignment):
    rng=np.random.default_rng(20261125)
    # All four categories are populated. Every category assignment can be
    # checked independently, with each tail choosing among three axes.
    data=base.chart(19,3,7)
    engine=Projections(data,3,1,1,assignment)
    checked=0
    for trial in range(30):
        scores=rng.normal(size=(6,4))
        states=engine.assign_scores(scores)
        got=sum((1-2*(s%2))*scores[i,s//2] for i,s in enumerate(states[1:]))
        best=-np.inf
        for labels in itertools.product(range(4),repeat=6):
            if [labels.count(c) for c in range(4)]!=engine.counts.tolist():continue
            score=sum(scores[i,0] if c==0 else -scores[i,0] if c==1 else max(scores[i,1:]) if c==2 else -min(scores[i,1:]) for i,c in enumerate(labels))
            best=max(best,score)
        if abs(got-best)>1e-9:raise ValueError('assignment differs from complete small enumeration')
        checked+=1
    # A separate resource DP allows all head counts while enforcing both signed
    # sums. Compare against every category word, not another assignment solver.
    free_engine=Projections(data,None,1,1,assignment)
    free_checked=0
    for trial in range(10):
        scores=rng.normal(size=(6,4))
        states=free_engine.assign_scores(scores)
        got=sum((1-2*(s%2))*scores[i,s//2] for i,s in enumerate(states[1:]))
        best=-np.inf
        for labels in itertools.product(range(4),repeat=6):
            if labels.count(0)-labels.count(1)!=1 or labels.count(2)-labels.count(3)!=1:continue
            score=sum(scores[i,0] if c==0 else -scores[i,0] if c==1 else max(scores[i,1:]) if c==2 else -min(scores[i,1:]) for i,c in enumerate(labels))
            best=max(best,score)
        if abs(got-best)>1e-9:raise ValueError('free-head projection differs from complete enumeration')
        if not np.array_equal(base.expand(states,data).sum(axis=1),free_engine.target):raise ValueError('free-head projection changed DC sums')
        free_checked+=1
    # Test covariance, conjugate symmetry and the Fourier normalization against
    # independent direct finite sums, rather than assuming NumPy's convention.
    data=base.chart(31,3,5);engine=Projections(data,1,-1,3,assignment)
    states,b=engine.discrete(rng.normal(size=(4,31)))
    if not np.array_equal(b.sum(axis=1),engine.target):raise ValueError('control projection changed counts')
    spectral=engine.spectral(np.fft.fft(b,axis=1))
    if np.max(np.abs(np.sum(np.abs(spectral)**2,axis=0)-31))>1e-9:raise ValueError('Fourier norm convention failed')
    values=np.fft.ifft(spectral,axis=1)
    if np.max(np.abs(values.imag))>1e-10:raise ValueError('conjugate symmetry failed')
    for x in range(31):
        if np.max(np.abs(values[:,5*x%31]-values[[0,3,1,2],x]))>1e-9:raise ValueError('three-cycle covariance failed')
    direct=np.array([sum(b[a,x]*np.exp(-2j*np.pi*k*x/31) for x in range(31)) for a in range(4) for k in range(31)]).reshape(4,31)
    if np.max(np.abs(direct-np.fft.fft(b,axis=1)))>1e-9:raise ValueError('FFT differs from direct Fourier sum')
    return {'complete_small_assignment_comparisons':checked,'complete_free_head_comparisons':free_checked,'DC_counts_checked':True,
            'spectral_norm_realness_covariance_checked':True,'direct_Fourier_comparison_checked':True}


def search(data,heads,head_net,tail_net,assignment,seconds,seed,beta,output,restart_after=2000,initial=None,orientation="discrete-first"):
    if not 0<beta<2 or seconds<=0:raise ValueError('positive duration and 0<beta<2 required')
    if orientation not in ('discrete-first','spectral-first'):raise ValueError('unknown projection orientation')
    engine=Projections(data,heads,head_net,tail_net,assignment)
    rng=np.random.default_rng(seed);start=time.monotonic();best=None
    iterations=0;restarts=0;stale=0
    if initial is None:
        _states,b=engine.discrete(rng.normal(size=(4,data['p'])))
    else:
        b=base.expand(initial,data)
        if not np.array_equal(b.sum(axis=1),engine.target):raise ValueError('initial sign counts differ')
    iterate=np.fft.fft(b,axis=1)
    while time.monotonic()-start<seconds:
        first=engine.spectral(iterate) if orientation=='spectral-first' else None
        reflected=2*first-iterate if first is not None else iterate
        states,b=engine.discrete(np.fft.ifft(reflected,axis=1).real)
        transformed=np.fft.fft(b,axis=1)
        numeric=np.rint(np.fft.ifft(np.sum(np.abs(transformed)**2,axis=0)).real).astype(np.int64)
        energy=int(numeric[data['shifts']]@numeric[data['shifts']])
        if best is None or energy<best['best_energy']:
            candidate=verify(states,data,engine.target)
            if candidate['best_energy']!=energy:raise ValueError('floating search objective differs from exact integers')
            best=candidate;stale=0
            best.update(p=data['p'],m=data['m'],mu=data['mu'],head_orbits=heads,
                        head_net=head_net,tail_net=tail_net,method='Douglas-Rachford',beta=beta,
                        head_stratum='all' if heads is None else heads,
                        occupied_head_orbits=sum(state//2==0 for state in states[1:]),
                        orientation=orientation,iterations=iterations,restarts=restarts,elapsed_seconds=time.monotonic()-start,
                        floating_iterates_are_only_search_guidance=True,completed_new_order_claimed=False)
            output.write_text(json.dumps(best,indent=2)+'\n')
            print(json.dumps({'best_energy':energy,'iterations':iterations,'seconds':best['elapsed_seconds']}),flush=True)
            if best['exact_T_design']:break
        if orientation=='discrete-first':
            projected=engine.spectral(2*transformed-iterate)
            iterate+=beta*(projected-transformed)
        else:
            iterate+=beta*(transformed-first)
        iterations+=1;stale+=1
        if stale>=restart_after:
            _states,b=engine.discrete(rng.normal(size=(4,data['p'])))
            iterate=np.fft.fft(b,axis=1);stale=0;restarts+=1
    best.update(total_iterations=iterations,total_restarts=restarts,total_seconds=time.monotonic()-start)
    output.write_text(json.dumps(best,indent=2)+'\n')
    return best


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--p',type=int,default=223);p.add_argument('--m',type=int,default=3);p.add_argument('--mu',type=int,default=183)
    p.add_argument('--heads',type=int,default=5);p.add_argument('--head-net',type=int,default=-5);p.add_argument('--tail-net',type=int,default=3)
    p.add_argument('--seconds',type=float,default=30);p.add_argument('--seed',type=int,default=20261125223)
    p.add_argument('--free-heads',action='store_true')
    p.add_argument('--orientation',choices=('discrete-first','spectral-first'),default='discrete-first')
    p.add_argument('--beta',type=float,default=.8);p.add_argument('--restart-after',type=int,default=2000)
    p.add_argument('--report',type=Path,required=True);p.add_argument('--workdir',type=Path);p.add_argument('--initial',type=Path)
    a=p.parse_args();folder=a.workdir or Path(tempfile.mkdtemp(prefix='t-projection-'))
    assignment,library=build_assignment(folder)
    control=controls(assignment);print(json.dumps({'controls':control}),flush=True)
    data=base.chart(a.p,a.m,a.mu);a.report.parent.mkdir(parents=True,exist_ok=True)
    initial=json.loads(a.initial.read_text())['states'] if a.initial else None
    report=search(data,None if a.free_heads else a.heads,a.head_net,a.tail_net,assignment,a.seconds,a.seed,a.beta,a.report,a.restart_after,initial,a.orientation)
    report.update(controls=control,assignment_source_sha256=hashlib.sha256((folder/'assignment-source.cpp').read_bytes()).hexdigest(),
                  runner_source_sha256=hashlib.sha256((folder/'projection-source.py').read_bytes()).hexdigest())
    a.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'exact_T_design':report['exact_T_design'],'best_energy':report['best_energy'],'iterations':report['total_iterations']}))

if __name__=='__main__':main()
